#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    QueryFormula, QueryMol, querymol,
    tidyformula, findformula,
    removehydrogens, inferaromaticity,
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
    # :recursive
    if a.key === :recursive && b.key === :recursive
        a == b && return true
        return hassubstructmatch(
            smartstomol(a.value), smartstomol(b.value),
            mandatory=Dict(1 => 1))
    elseif a.key === :recursive && !(b.key in (:not, :or))
        amol = smartstomol(a.value)
        return issubset(nodeattr(amol, 1).query, b)
    elseif b.key === :recursive
        return false
    end
    if !(a.key in (:and, :or) || b.key in (:and, :or))
        # :not
        a.key === :not && b.key === :not && return a == b
        a.key === :not && return false
        b.key === :not && a.key !== :recursive && return a.key == b.value.key && a.value != b.value.value
        # descriptors
        return a == b
    end
    # :and, :or
    aset = a.key in (:and, :or) ? a.value : [a]
    bset = b.key in (:and, :or) ? b.value : [b]
    akeys = Set(e.key for e in aset)
    bkeys = Set(e.key for e in bset)
    if a.key === :and
        if b.key === :and && length(bkeys) == 1 && collect(bkeys)[1] === :not
            return any(issubset(fml, b) for fml in a.value)
        else
            amap = Dict(i => v for (i, v) in enumerate(aset))
            bmap = Dict(i => v for (i, v) in enumerate(b.key === :and ? bset : [b]))
            func = (x, y) -> issubset(amap[x], bmap[y])
            return maxcard(keys(amap), keys(bmap), func) == length(bmap)
        end
    elseif b.key === :or
        amap = Dict(i => v for (i, v) in enumerate(a.key === :or ? aset : [a]))
        bmap = Dict(i => v for (i, v) in enumerate(bset))
        func = (x, y) -> issubset(amap[x], bmap[y])
        return maxcard(keys(amap), keys(bmap), func) == length(amap)
    elseif b.key === :and && length(bkeys) == 1 && collect(bkeys)[1] === :not
        for bfml in bset
            for afml in aset
                issubset(afml, bfml) || return false
            end
        end
        return true
    elseif a.key === :or && b.key === :not
        if length(akeys) == 1 && collect(akeys)[1] == b.value.key
            amap = Dict(i => v for (i, v) in enumerate(aset))
            bmap = Dict(i => v for (i, v) in enumerate(bset))
            func = (x, y) -> issubset(amap[x], bmap[y])
            return maxcard(keys(amap), keys(bmap), func) == 1
        end
    end
    return false
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
        if fml.key === :and && length(mono) == 1 && !(collect(mono)[1] in (:not, :recursive))
            # SMARTS primitives are disjoint, so the intersection should be an empty set.
            # NOTE: except for ! and $() queries
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
    removehydrogens(mol::QueryMol) -> QueryMol

Return the molecular query with hydrogen nodes removed.
"""
function removehydrogens(qmol::QueryMol)
    # count H nodes and mark H nodes to remove
    hnodes = Set{Int}()
    hcntarr = zeros(Int, nodecount(qmol))
    for n in 1:nodecount(qmol)
        nq = nodeattr(qmol, n).query
        issubset(nq, QueryFormula(:atomsymbol, :H)) || continue
        degree(qmol, n) == 1 || throw(ErrorException("Invalid hydrogen valence"))
        adj = iterate(adjacencies(qmol, n))[1]
        hcntarr[adj] += 1
        push!(hnodes, n)
    end
    qmol_ = deepcopy(qmol)
    heavynodes = setdiff(nodeset(qmol), hnodes)
    for n in heavynodes
        nq = nodeattr(qmol, n).query
        # [!#1] => *
        noth = QueryFormula(:not, QueryFormula(:atomsymbol, :H))
        if nq == noth
            newq = QueryFormula(:any, true)
        elseif nq.key === :and && noth in nq.value
            if length(nq.value) == 2
                newq = collect(setdiff(nq.value, [noth]))[1]
            else
                newq = QueryFormula(:and, setdiff(nq.value, [noth]))
            end
        else
            newq = nq
        end
        # consider H nodes
        adjhnfmls = collect(0:(hcntarr[n] - 1))
        if !isempty(adjhnfmls)
            nfmls = [QueryFormula(:not, QueryFormula(:hydrogenconnected, i)) for i in adjhnfmls]
            newq = QueryFormula(:and, Set([newq, nfmls...]))
        end
        setnodeattr!(qmol_, n, SmartsAtom(tidyformula(newq)))
    end
    return querymol(nodesubgraph(qmol_, heavynodes))
end


"""
    inferatomaromaticity(qmol::QueryMol)

Infer aromaticity of atoms and bonds, then return more specific query in the aspect of aromaticity.
"""
function inferaromaticity(qmol::QueryMol)
    qmol_ = deepcopy(qmol)
    for n in 1:nodecount(qmol)
        nq = nodeattr(qmol, n).query
        issubset(nq, QueryFormula(:isaromatic, true)) && continue
        issubset(nq, QueryFormula(:isaromatic, false)) && continue
        # by topology query (!R, !r)
        if issubset(nq, QueryFormula(:sssrcount, 0))
            newq = QueryFormula(:and, Set([nq, QueryFormula(:isaromatic, false)]))
            setnodeattr!(qmol_, n, SmartsAtom(tidyformula(newq)))
            continue
        end
        # by atom symbol
        canbearom = [:B, :C, :N, :O, :P, :S, :As, :Se]
        notaromfml = QueryFormula(:and, Set([
            QueryFormula(:not, QueryFormula(:atomsymbol, a)) for a in canbearom]))
        if issubset(nq, notaromfml)
            newq = QueryFormula(:and, Set([nq, QueryFormula(:isaromatic, false)]))
            setnodeattr!(qmol_, n, SmartsAtom(tidyformula(newq)))
            continue
        end
        # by explicitly non-/aromatic incidences
        noincacc = QueryFormula(:and, Set([
            QueryFormula(:not, QueryFormula(:atomsymbol, a)) for a in [:C, :N, :B]]))
        minacc = issubset(nq, noincacc) ? 0 : 1
        nonaromcnt = 0
        # hydrogen query
        if issubset(nq, QueryFormula(:and, Set([
            QueryFormula(:not, QueryFormula(:hydrogenconnected, 0)),
            QueryFormula(:not, QueryFormula(:hydrogenconnected, 1))
        ])))
            nonaromcnt += 2  # 2 is enough to be nonaromcnt > minacc
        elseif issubset(nq, QueryFormula(:not, QueryFormula(:hydrogenconnected, 0)))
            nonaromcnt += 1
        end
        # incidences
        hasarombond = false
        for inc in incidences(qmol, n)
            eq = edgeattr(qmol_, inc).query
            if issubset(eq, QueryFormula(:isaromaticbond, true))
                hasarombond = true
                break
            end
            if issubset(eq, QueryFormula(:bondorder, 2))
                # C=O special case
                adjq = nodeattr(qmol, neighbors(qmol, n)[inc]).query
                if issubset(adjq, QueryFormula(:atomsymbol, :O))
                    continue
                end
            end
            if (issubset(eq, QueryFormula(:isaromaticbond, false))
                    || issubset(eq, QueryFormula(:isringbond, false)))
                if issubset(eq, QueryFormula(:not, QueryFormula(:bondorder, 1)))
                    nonaromcnt += 2
                else
                    nonaromcnt += 1
                end
            end
        end
        if hasarombond
            newq = QueryFormula(:and, Set([nq, QueryFormula(:isaromatic, true)]))
            setnodeattr!(qmol_, n, SmartsAtom(tidyformula(newq)))
        elseif nonaromcnt > minacc
            newq = QueryFormula(:and, Set([nq, QueryFormula(:isaromatic, false)]))
            setnodeattr!(qmol_, n, SmartsAtom(tidyformula(newq)))
        end
    end
    # by Huckel rule
    aromf = QueryFormula(:isaromaticbond, true)
    aors = QueryFormula(:or, Set([
        aromf,
        QueryFormula(:and, Set([
            QueryFormula(:bondorder, 1),
            QueryFormula(:isaromaticbond, false)
        ]))
    ]))
    aord = QueryFormula(:or, Set([
        aromf,
        QueryFormula(:and, Set([
            QueryFormula(:bondorder, 2),
            QueryFormula(:isaromaticbond, false)
        ]))
    ]))
    for ring in sssr(qmol)
        ringedges = edgeset(nodesubgraph(qmol, ring))
        pcnt = 0
        for n in ring
            nq = nodeattr(qmol, n).query
            if issubset(nq, QueryFormula(:isaromatic, true))
                pcnt += 1
                continue
            end
            rincs = collect(intersect(incidences(qmol, n), ringedges))
            uq = edgeattr(qmol, rincs[1]).query
            vq = edgeattr(qmol, rincs[2]).query
            if uq == aors && vq == aord || vq == aors && uq == aord
                if issubset(nq, QueryFormula(:or, Set([QueryFormula(:atomsymbol, a) for a in [:B, :C, :N, :P, :As]])))
                    pcnt += 1
                    continue
                end
            elseif uq == aors && vq == aors
                if issubset(nq, QueryFormula(:or, Set([QueryFormula(:atomsymbol, a) for a in [:N, :O, :P, :S, :As, :Se]])))
                    pcnt += 2
                    continue
                elseif issubset(nq, QueryFormula(:atomsymbol, :C))
                    outer = collect(setdiff(incidences(qmol, n), ringedges))
                    if length(outer) == 1
                        outerq = edgeattr(qmol, outer[1]).query
                        oadjq = nodeattr(qmol, neighbors(qmol, n)[outer[1]]).query
                        if issubset(outerq, QueryFormula(:bondorder, 2)) && issubset(oadjq, QueryFormula(:atomsymbol, :O))
                            continue
                        end
                    end
                end
            end
            pcnt = 0
            break
        end
        if pcnt % 4 == 2
            for n in ring
                nq = nodeattr(qmol, n).query
                newq = QueryFormula(:and, Set([nq, QueryFormula(:isaromatic, true)]))
                setnodeattr!(qmol_, n, SmartsAtom(tidyformula(newq)))
            end
            for e in ringedges
                setedgeattr!(qmol_, e, SmartsBond(aromf))
            end
        end
    end
    return qmol_
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
        rcd["parsed"] = smartstomol(rcd["query"])
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
function filter_queries(qr::DictDiGraph, mol::GraphMol; filtering=true)
    matched = Set{Int}()
    qrc = deepcopy(qr)
    for n in reversetopologicalsort(qrc)
        rcd = nodeattr(qrc, n)
        if filtering
            if !issubset(successors(qrc, n), matched)  # query containment filter
                continue
            end
        end
        matches = collect(substructmatches(mol, rcd["parsed"]))
        if !isempty(matches)
            push!(matched, n)
            rcd["matched"] = Set([sort(collect(keys(m))) for m in matches])
        end
    end
    return dictdigraph(nodesubgraph(qrc, matched))
end