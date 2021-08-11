#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    QueryFormula, QueryMol, querymol,
    tidyformula, findformula,
    inferatomaromaticity!, removehydrogens


struct QueryFormula
    key::Symbol
    value::Any
end


Base.:(==)(a::QueryFormula, b::QueryFormula) = a.key == b.key && a.value == b.value

function Base.in(a::QueryFormula, b::Set{QueryFormula})
    # Necessary for set operations
    for i in b
        a == i && return true
    end
    return false
end


"""
    issubset(a::QueryFormula, b::QueryFormula) -> Bool

Check if fml1 contains fml2 (that is, all the query results of fml1 is included in the results of fml2)
"""
function Base.issubset(a::QueryFormula, b::QueryFormula)
    b == QueryFormula(:any, true) && return true
    if a.key === :and
        amap = Dict(i => v for (i, v) in enumerate(a.value))
        bset = b.key === :and ? b.value : [b]
        bmap = Dict(i => v for (i, v) in enumerate(bset))
        func = (x, y) -> issubset(amap[x], bmap[y])
        return maxcard(keys(amap), keys(bmap), func) == length(bmap)
    elseif b.key === :or
        bmap = Dict(i => v for (i, v) in enumerate(b.value))
        aset = a.key === :or ? a.value : [a]
        amap = Dict(i => v for (i, v) in enumerate(aset))
        func = (x, y) -> issubset(amap[x], bmap[y])
        return maxcard(keys(amap), keys(bmap), func) == length(amap)
    else
        return a == b
    end
end



struct QueryMol{A<:QueryAtom,B<:QueryBond} <: OrderedGraph
    neighbormap::Vector{Dict{Int,Int}}
    edges::Vector{Tuple{Int,Int}}
    nodeattrs::Vector{A}
    edgeattrs::Vector{B}
    cache::Dict{Symbol,Any}
    attributes::Dict{Symbol,Any}
    connectivity::Vector{Vector{Int}}
end

"""
    querymol() -> QueryMol

Generate empty `QueryMol`.
"""
querymol(
    ::Type{A}, ::Type{B}
) where {A<:QueryAtom,B<:QueryBond} = QueryMol{A,B}(
    [], [], [], [], Dict(), Dict(), [])

"""
    querymol(atoms::Vector{Atom}, bonds::Vector{Bond}) -> GraphMol

Generate `QueryMol` that has the given atom objects and edge objects.
"""
function querymol(edges, atoms::Vector{A}, bonds::Vector{B},
        connectivity::Vector{Vector{Int}}) where {A<:QueryAtom,B<:QueryBond}
    nbrmap = [Dict{Int,Int}() for i in 1:length(atoms)]
    edges = collect(edges)
    for (i, (u, v)) in enumerate(edges)
        nbrmap[u][i] = v
        nbrmap[v][i] = u
    end
    return QueryMol(
        nbrmap, edges, atoms, bonds,
        Dict{Symbol,Any}(), Dict{Symbol,Any}(), connectivity)
end


"""
    querymol(mol::SubgraphView{QueryMol}) -> QueryMol

Generate a new `QueryMol` from a substructure view.

Graph property caches and attributes are not inherited.
"""
function querymol(view::SubgraphView)
    newg = querymol(nodeattrtype(view), edgeattrtype(view))
    nkeys = sort(collect(nodeset(view)))
    ekeys = sort(collect(edgeset(view)))
    nmap = Dict{Int,Int}()
    for (i, n) in enumerate(nkeys)
        nmap[n] = i
        push!(newg.nodeattrs, nodeattr(view, n))
        push!(newg.neighbormap, Dict())
    end
    for (i, e) in enumerate(ekeys)
        (oldu, oldv) = getedge(view, e)
        u = nmap[oldu]
        v = nmap[oldv]
        push!(newg.edges, (u, v))
        push!(newg.edgeattrs, edgeattr(view, e))
        newg.neighbormap[u][i] = v
        newg.neighbormap[v][i] = u
    end
    return newg
end



"""
tidyformula(fml::QueryFormula) -> QueryFormula

Return juxtaposed formulae (ex. :and => (A, :and => (B, C)) -> :and => (A, B, C)).
"""
function tidyformula(fml::QueryFormula)
    fml.key in (:and, :or) || return fml
    # TODO: :not and :recursive formula
    childs = Set{QueryFormula}()
    # Association
    for child in fml.value
        if child.key === fml.key
            union!(childs, tidyformula(child).value)
        else
            union!(childs, [child])
        end
    end
    length(childs) == 1 && return collect(childs)[1]
    # Distribution
    ckey = fml.key === :and ? :or : :and
    cc = collect(childs)
    bin = cc[1].key == ckey ? copy(cc[1].value) : Set([cc[1]])
    for child in cc[2:end]
        intersect!(bin, child.key == ckey ? child.value : [child])
    end
    isempty(bin) && return QueryFormula(fml.key, childs)
    updated = Set{QueryFormula}()
    for child in childs
        cdiff = setdiff(child.key == ckey ? child.value : [child], bin)
        if isempty(cdiff)
            # Absorption
            return length(bin) == 1 ? collect(bin)[1] : QueryFormula(ckey, bin)
        elseif length(cdiff) == 1
            push!(updated, collect(cdiff)[1])
        else
            push!(updated, QueryFormula(ckey, cdiff))
        end
    end
    return tidyformula(QueryFormula(ckey, Set([
        bin...,
        QueryFormula(fml.key, updated)
    ])))
end


"""
    findformula(q::QueryFormula, key::Symbol) -> Any

Find value of the given key.

If `deep=true` is set, it will also search under `:or` sets to match multiple formulae and
return the value after applying the `aggregate` function to them.
"""
function findformula(q::QueryFormula, key::Symbol; deep=false, aggregate=minimum)
    if q.key == key
        return q.value
    elseif q.key === :and
        for child in q.value
            if child.key == key
                return child.value
            elseif deep && child.key === :or
                gc = findformula(child, key, deep=deep, aggregate=aggregate)
                gc === nothing || return gc
            end
        end
    elseif deep && q.key === :or
        values = []
        for child in q.value
            child.key == key || return
            push!(values, child.value)
        end
        return aggregate(values)
    end
    return
end


"""
    inferatomaromaticity!(qmol::QueryMol)

Infer aromaticity of atoms from connected bonds.

This directly add :isaromatic formulas to query atoms in the query molecule.
"""
function inferatomaromaticity!(qmol::QueryMol)
    for n in 1:nodecount(qmol)
        arom = false
        hcnt = 0
        heavyorder = 0
        for (inc, adj) in neighbors(qmol, n)
            iq = edgeattr(qmol, inc).query
            if findformula(iq, :isaromaticbond) === true
                arom = true
                break
            end
            aq = nodeattr(qmol, adj).query
            if findformula(aq, :atomsymbol) === :H
                hcnt += 1
            else
                heavyorder += findformula(iq, :bondorder, deep=true)
            end
        end
        nq = nodeattr(qmol,n).query
        if arom
            newq = QueryFormula(:and, Set([nq, QueryFormula(:isaromatic, true)]))
            setnodeattr!(qmol, n, SmartsAtom(tidyformula(newq)))
            continue
        end
        hqcnt = findformula(nq, :hydrogenconnected, deep=true)
        if hqcnt !== nothing
            hcnt = hqcnt  # ignore hydrogen nodes if H query is already set
        end
        acceptable = Dict(:O => 0, :N => 1, :C => 1, :S => 0, :P => 1, :B => 1)
        if heavyorder + hcnt > get(acceptable, findformula(nq, :atomsymbol), 0)
            newq = QueryFormula(:and, Set([nq, QueryFormula(:isaromatic, false)]))
            setnodeattr!(qmol, n, SmartsAtom(tidyformula(newq)))
        end
    end
end


"""
    removehydrogens(mol::QueryMol) -> QueryMol

Return the molecular query with hydrogen nodes removed.

Hydrogen nodes would be represented as `[*Hx]` in the connected heavy atoms.
(e.g. CN[H] -> C[N;H1,H2,H3,H4])
"""
function removehydrogens(mol::QueryMol)
    MAX_H_COUNT = 4
    hs = Set{Int}()
    hcntmap = Dict()
    for (i, atom) in enumerate(nodeattrs(mol))
        h = atom.query == QueryFormula(:atomsymbol, :H)
        h2 = atom.query == QueryFormula(
            :and, Set([QueryFormula(:atomsymbol, :H), QueryFormula(:isaromatic, false)]))
        (h || h2) || continue
        @assert degree(mol, i) == 1
        adj = iterate(adjacencies(mol, i))[1]
        if !haskey(hcntmap, adj)
            hcntmap[adj] = 0
        end
        hcntmap[adj] += 1
        push!(hs, i)
    end
    newmol = deepcopy(mol)
    for (i, hcnt) in hcntmap
        q = nodeattr(mol, i).query
        findformula(q, :hydrogenconnected, deep=true) === nothing || continue  # ignore if H is already set
        heavycnt = degree(mol, i) - hcnt
        maxh = MAX_H_COUNT - heavycnt
        if maxh <= hcnt
            hf = QueryFormula(:hydrogenconnected, hcnt)
        else
            hf = QueryFormula(
                :or, Set([QueryFormula(:hydrogenconnected, n) for n in hcnt:maxh]))
        end
        newq = QueryFormula(:and, Set([q, hf]))
        setnodeattr!(newmol, i, SmartsAtom(tidyformula(newq)))
    end
    ns = setdiff(nodeset(mol), hs)
    return querymol(nodesubgraph(newmol, ns))
end
