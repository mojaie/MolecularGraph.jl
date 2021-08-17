#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    QueryFormula, QueryMol, querymol,
    tidyformula, findformula, findallformula,
    convertnotquery, convertnotquery!,
    inferatomaromaticity!, removehydrogens,
    recursiveatommatch,
    query_relationship, filter_queries


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
    elseif a.key === :or && b.key === :not
        common = Set([c.key for c in a.value])
        length(common) == 1 || return false
        collect(common)[1] == b.value.key || return false
        return !(b.value.value in Set([c.value for c in a.value]))
    elseif a.key === :recursive
        if b.key === :recursive
            return hassubstructmatch(
                smartstomol(a.value), smartstomol(b.value),
                mandatory=Dict(1 => 1)
            )
        else
            amol = smartstomol(a.value)
            return issubset(nodeattr(amol, 1).query, b)
        end
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

Return tidy formulae.

- associative formulae will be juxtaposed
  (ex. :and => (A, :and => (B, C)) -> :and => (A, B, C))
- distributive formulae will be factored out
  (ex. :or => (:and => (A, B), :and => (A, C)) -> :and => (A, :or => (B, C)))
- Absorption
  (ex. :and => (A, :or => (A, B)) -> A
- `:any` absorbs everything
  (ex. :and => (:any => true, A) -> A, :or => (:any => true, A) -> :any => true
- `:not` would be inverted if possible
  (ex. :not => (:A => true) -> :A => false)
- `:not` has the highest precedence in SMARTS, but only in the case like [!C],
  De Morgan's law will be applied to remove `:and` under `:not`.
  (ex. :not => (:and => (:atomsymbol => :C, :isaromatic => false)
   -> :or => (:not => (:atomsymbol => :C), isaromatic => true)
"""
function tidyformula(fml::QueryFormula)
    # not
    if fml.key === :not
        child = fml.value
        if child.key === :and  # only the cases like [!C]
            return QueryFormula(:or, Set([
                tidyformula(QueryFormula(:not, c)) for c in child.value
            ]))
        elseif typeof(child.value) === Bool
            return QueryFormula(child.key, !child.value)
        end
    end
    fml.key in (:and, :or) || return fml
    childs = Set{QueryFormula}()
    # Association
    for child in fml.value
        cfml = tidyformula(child)
        if cfml.key === :any
            # Absorption
            (fml.key === :and) == cfml.value && continue
            return QueryFormula(:any, cfml.value)
        elseif cfml.key === fml.key
            union!(childs, cfml.value)
        else
            union!(childs, [cfml])
        end
    end
    length(childs) == 1 && return collect(childs)[1]
    # Distribution
    ckey = fml.key === :and ? :or : :and
    cc = collect(childs)
    bin = cc[1].key == ckey ? copy(cc[1].value) : Set([cc[1]])
    mono = Set([b.key for b in bin])
    for child in cc[2:end]
        @assert child.key !== fml.key  # already associated
        elems = child.key === ckey ? child.value : [child]
        intersect!(bin, elems)
        union!(mono, Set([e.key for e in elems]))
    end
    if isempty(bin)
        if fml.key === :and && length(mono) == 1
            # SMARTS primitives are disjoint, so the intersection should be an empty set.
            return QueryFormula(:any, false)
        else
            return QueryFormula(fml.key, childs)
        end
    end
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
    return tidyformula(QueryFormula(ckey, Set([bin..., QueryFormula(fml.key, updated)])))
end


"""
    findformula(q::QueryFormula, key::Symbol) -> Any

Find value of the given key.

If `deep=true` is set, it will also search under `:or` sets to match multiple formulae and
return the value after applying the `aggregate` function to them.
"""
function findformula(q::QueryFormula, key::Symbol; deep=false, aggregate=minimum)
    if q.key == key
        return aggregate([q.value])
    elseif q.key === :and
        for child in q.value
            if child.key == key
                return aggregate([child.value])
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
    findallformula(q::QueryFormula, key::Symbol) -> Any

Find all formula values that matched the given key regardless of logical operator roles.
"""
function findallformula(q::QueryFormula, key::Symbol)
    if q.key == key
        return [q.value]
    elseif q.key === :not
        return findallformula(q.value, key)
    elseif q.key in (:and, :or)
        values = []
        for child in q.value
            append!(values, findallformula(child, key))
        end
        return values
    end
    return []
end



NOTCONV = Dict(
    :isaromatic => Set([true, false]),
    :hydrogenconnected => Set([0, 1, 2, 3, 4]),
    :sssrcount => Set([0, 1, 2, 3]),
    :isringbond => Set([true, false]),
    :isaromaticbond => Set([true, false]),
    :stereo => Set([:up, :down, :unspecified])
)

"""
    convertnotquery(mol::QueryMol) -> QueryMol

Convert :not query to :or if possible.

- [!#1] would be *. Typical substructure mining tasks do not care stereochemistry,
so all hydrogen nodes can be removed.)
- following `tidyformula` is necessary.
"""
function convertnotquery(fml::QueryFormula)
    if fml.key === :not
        cfml = fml.value
        cfml == QueryFormula(:atomsymbol, :H) && return QueryFormula(:any, true)
        if haskey(NOTCONV, cfml.key)
            rest = setdiff(NOTCONV[cfml.key], [cfml.value])
            isempty(rest) && return QueryFormula(:any, false)
            length(rest) == 1 && return QueryFormula(cfml.key, collect(rest)[1])
            return QueryFormula(:or, Set([QueryFormula(cfml.key, v) for v in rest]))
        else
            return fml
        end
    end
    fml.key in (:and, :or) || return fml
    childs = Set{QueryFormula}(convertnotquery.(fml.value))
    return QueryFormula(fml.key, childs)
end

function convertnotquery!(qmol::QueryMol)
    nodes = []
    for n in 1:nodecount(qmol)
        na = nodeattr(qmol, n).query
        push!(nodes, SmartsAtom(tidyformula(convertnotquery(na))))
    end
    edges = []
    for e in 1:edgecount(qmol)
        ea = edgeattr(qmol, e).query
        push!(edges, SmartsBond(tidyformula(convertnotquery(ea))))
    end
    empty!(qmol.nodeattrs)
    append!(qmol.nodeattrs, nodes)
    empty!(qmol.edgeattrs)
    append!(qmol.edgeattrs, edges)
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
                heavyorder += something(findformula(iq, :bondorder, deep=true), 1)
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

- Hydrogen nodes would be represented as `[*Hx]` in the connected heavy atoms.
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


function recursiveatommatch(qmol1::QueryMol, qmol2::QueryMol)
    return function (a1, a2)
        a1q = nodeattr(qmol1, a1).query
        a2q = nodeattr(qmol2, a2).query
        issubset(a1q, a2q) && return true
        recfmls= findformula(a2q, :recursive, deep=true, aggregate=collect)
        recfmls === nothing && return false
        for r in recfmls
            hassubstructmatch(
                qmol1, smartstomol(r), mandatory=Dict(a1 => 1)) && return true
        end
        return false
    end
end



const DEFAULT_QUERY_RELATIONS = let
    qrfile = joinpath(dirname(@__FILE__), "../../assets/const/default_query_relations.yaml")
    include_dependency(qrfile)
    qrfile
end



struct DictDiGraph <: OrderedDiGraph
    # TODO: should be moved to Graph module
    outneighbormap::Vector{Dict{Int,Int}}
    inneighbormap::Vector{Dict{Int,Int}}
    edges::Vector{Tuple{Int,Int}}
    nodeattrs::Vector{Dict}
    edgeattrs::Vector{Dict}
end


"""
    dictdigraph(view::DiSubgraphView{DictDiGraph}) -> DictDiGraph

Generate a new `DictDiGraph` from a substructure view.

Graph property caches and attributes are not inherited.
"""
function dictdigraph(view::DiSubgraphView{DictDiGraph})
    newg = DictDiGraph([], [], [], [], [])
    nkeys = sort(collect(nodeset(view)))
    ekeys = sort(collect(edgeset(view)))
    nmap = Dict{Int,Int}()
    for (i, n) in enumerate(nkeys)
        nmap[n] = i
        push!(newg.nodeattrs, nodeattr(view, n))
        push!(newg.outneighbormap, Dict())
        push!(newg.inneighbormap, Dict())
    end
    for (i, e) in enumerate(ekeys)
        (oldu, oldv) = getedge(view, e)
        u = nmap[oldu]
        v = nmap[oldv]
        push!(newg.edges, (u, v))
        push!(newg.edgeattrs, edgeattr(view, e))
        newg.outneighbormap[u][i] = v
        newg.inneighbormap[v][i] = u
    end
    return newg
end


"""
    query_relationship(;sourcefile=DEFAULT_QUERY_RELATIONS) -> DictDiGraph

Generate query relationship diagram.
"""
function query_relationship(;sourcefile=DEFAULT_QUERY_RELATIONS)
    graph = DictDiGraph([], [], [], [], [])
    keys = Dict()
    for (i, rcd) in enumerate(YAML.load(open(sourcefile)))
        addnode!(graph, rcd)
        keys[rcd["key"]] = i
    end
    for rcd in nodeattrs(graph)
        if haskey(rcd, "isa")
            for e in rcd["isa"]
                addedge!(graph, keys[rcd["key"]], keys[e], Dict("relation" => "isa"))
            end
        end
        if haskey(rcd, "has")
            for e in rcd["has"]
                addedge!(graph, keys[rcd["key"]], keys[e], Dict("relation" => "has"))
            end
        end
    end
    return graph
end


"""
    filter_queries(qr::DictDiGraph, mol::GraphMol) -> DictDiGraph

Filter query relationship diagram by the given molecule.
The filtered diagram represents query relationship that the molecule have.
"""
function filter_queries(qr::DictDiGraph, mol::GraphMol)
    matched = Set{Int}()
    qrc = deepcopy(qr)
    for n in reversetopologicalsort(qrc)
        rcd = nodeattr(qrc, n)
        if !issubset(successors(qrc, n), matched)
            continue
        end
        matches = collect(substructmatches(mol, smartstomol(rcd["query"])))
        if !isempty(matches)
            push!(matched, n)
            rcd["matched"] = Set([sort(collect(keys(m))) for m in matches])
        end
    end
    return dictdigraph(nodesubgraph(qrc, matched))
end