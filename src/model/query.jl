#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    QueryLiteral, QueryTruthTable,
    statement, any_query, not_query, and_query, or_query
    # removehydrogens, inferaromaticity,
    # query_relationship, filter_queries


struct QueryLiteral
    operator::Symbol  # :eq, :gt?, :lt? ...
    key::Symbol
    value::Union{Symbol,Int,String,Bool,Nothing}
end

QueryLiteral(key) = QueryLiteral(:eq, key, true)
QueryLiteral(key, value) = QueryLiteral(:eq, key, value)

Base.isless(q::QueryLiteral, r::QueryLiteral
    ) = q.key < r.key || q.operator < r.operator || string(q.value) < string(r.value)


struct QueryTruthTable
    table::Matrix{Union{Bool,Nothing}}
    props::Vector{QueryLiteral}
end


"""
    QueryTruthTable(fml::Function, props::Vector{QueryLiteral}) -> QueryTruthTable

Construct QueryTruthTable from boolian function and QueryLiteral vector.

Maybe faster than query generator functions (e.g. and_query, or_query).  
- props: QueryLiteral vector. The properties should be unique.
- function: function that takes a vector whose size is `length(props)`
  that corresponds to each property variables and returns true or false.
"""
function QueryTruthTable(fml::Function, props::Vector{QueryLiteral})
    nlit = length(props)
    table = iszero.(transpose(hcat([digits(i, base=2, pad=nlit)[1:nlit] for i in 1:(2^nlit)]...)))
    values = Vector{Union{Bool,Nothing}}[]
    for i in 1:size(table, 1)
        vals = table[i, :]
        res = fml(vals)
        if res == 1
            push!(vals, res)
            push!(values, vals)
        end
    end
    return QueryTruthTable(transpose(hcat(values...)), props)
end

QueryTruthTable(fml::Function, props::Vector{T}
    ) where T <: Tuple = QueryTruthTable(fml, [QueryLiteral(p...) for p in props])


function issubset_arr(q, r)
    matched = true
    for i in 1:length(q)
        q[i] == r[i] && continue
        isnothing(q[i]) && continue
        isnothing(r[i]) && continue
        matched = false
        break
    end
    matched
end

function Base.issubset(q::QueryTruthTable, r::QueryTruthTable)
    # tautology
    """
    q \\ r  T other F
    T      T   F   F
    other  T   ?   F
    F      T   T   T
    """
    isempty(q.props) && (q.table[1] || return true)
    isempty(r.props) && (r.table[1] && return true)
    (isempty(q.props) || isempty(r.props)) && return false
    # TODO: more efficient implementation
    unn = sort(collect(union(Set(q.props), Set(r.props))))
    qtbl = Matrix{Union{Bool,Nothing}}(nothing, size(q.table, 1), length(unn) + 1)
    rtbl = Matrix{Union{Bool,Nothing}}(nothing, size(r.table, 1), length(unn) + 1)
    for (i, k) in enumerate(unn)
        for (m, p) in enumerate(q.props)
            if k == p  # TODO: recursive and disjoint check
                qtbl[:, i] = q.table[:, m]
                break
            end
        end
        qtbl[:, end] = q.table[:, end]
        for (m, p) in enumerate(r.props)
            if k == p
                rtbl[:, i] = r.table[:, m]
                break
            end
        end
        rtbl[:, end] = r.table[:, end]
    end
    for i in 1:size(qtbl, 1)
        matched = 0
        nreq = 2^sum(isnothing.(qtbl[i, :]))
        for j in 1:size(rtbl, 1)
            if issubset_arr(qtbl[i, :], rtbl[j, :])
                matched += 1
                matched == nreq && break
            end
        end
        matched == nreq || return false
    end
    return true
end

# TODO: expensive. just for testing
Base.:(==)(q::QueryTruthTable, r::QueryTruthTable) = issubset(q, r) && issubset(r, q)

# Query generation
statement(key, value=true, operator=:eq) = (v -> v[1], [(operator, key, value)])
any_query(b::Bool) = (v -> b, Tuple{Symbol}[])
not_query(q) = (v -> ~q[1](v), q[2])
function and_query(qs; cond=all)
    funcs = Function[]
    props = Union{Tuple{Symbol},Tuple{Symbol,Any},Tuple{Symbol,Symbol,Any}}[]
    vpos = Vector{Int}[]  # positions of vars passed to downstream
    for q in qs
        func, ps = q
        push!(funcs, func)
        idxs = Int[]
        for p in ps
            if p in props  # props should be unique
                idx = findfirst(x -> x == p, props)
                push!(idxs, idx)
            else
                push!(props, p)
                push!(idxs, length(props))
            end
        end
        push!(vpos, idxs)
    end
    func = v -> cond(f(v[vpos[i]]) for (i, f) in enumerate(funcs))
    return (func, props)
end
or_query(qs) = and_query(qs, cond=any)



"""
    removehydrogens(mol::QueryMol) -> QueryMol

Return the molecular query with hydrogen nodes removed.

function removehydrogens(qmol::QueryMol)
    # count H nodes and mark H nodes to remove
    hnodes = Set{Int}()
    hcntarr = zeros(Int, nodecount(qmol))
    for n in 1:nodecount(qmol)
        nq = nodeattr(qmol, n).query
        issubset(nq, QueryFormula(:symbol, :H), eval_recursive=false) || continue
        degree(qmol, n) == 1 || throw(ErrorException("Invalid hydrogen valence"))
        adj = iterate(adjacencies(qmol, n))[1]
        hcntarr[adj] += 1
        push!(hnodes, n)
    end
    qmol_ = deepcopy(qmol)
    heavynodes = setdiff(nodeset(qmol), hnodes)
    for n in heavynodes
        nq = nodeattr(qmol, n).query
        # no longer H nodes exist, so [!#1] would be ignored
        noth = QueryFormula(:not, QueryFormula(:symbol, :H))
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
        # consider H nodes as a :hydrogenconnected query
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

"""
    inferatomaromaticity(qmol::QueryMol)

Infer aromaticity of atoms and bonds, then return more specific query in the aspect of aromaticity.

function inferaromaticity(qmol::QueryMol)
    qmol_ = deepcopy(qmol)
    for n in 1:nodecount(qmol)
        nq = nodeattr(qmol, n).query
        issubset(nq, QueryFormula(:isaromatic, true), eval_recursive=false) && continue
        issubset(nq, QueryFormula(:isaromatic, false), eval_recursive=false) && continue
        # by topology query (!R, !r)
        if issubset(nq, QueryFormula(:sssrcount, 0), eval_recursive=false)
            newq = QueryFormula(:and, Set([nq, QueryFormula(:isaromatic, false)]))
            setnodeattr!(qmol_, n, SmartsAtom(tidyformula(newq)))
            continue
        end
        # by atom symbol
        canbearom = [:B, :C, :N, :O, :P, :S, :As, :Se]
        notaromfml = QueryFormula(:and, Set([
            QueryFormula(:not, QueryFormula(:symbol, a)) for a in canbearom]))
        if issubset(nq, notaromfml, eval_recursive=false)
            newq = QueryFormula(:and, Set([nq, QueryFormula(:isaromatic, false)]))
            setnodeattr!(qmol_, n, SmartsAtom(tidyformula(newq)))
            continue
        end
        # by explicitly non-/aromatic incidences
        noincacc = QueryFormula(:and, Set([
            QueryFormula(:not, QueryFormula(:symbol, a)) for a in [:C, :N, :B]]))
        minacc = issubset(nq, noincacc, eval_recursive=false) ? 0 : 1
        nonaromcnt = 0
        # hydrogen query
        if issubset(nq, QueryFormula(:and, Set([
            QueryFormula(:not, QueryFormula(:hydrogenconnected, 0)),
            QueryFormula(:not, QueryFormula(:hydrogenconnected, 1))
        ])), eval_recursive=false)
            nonaromcnt += 2  # 2 is enough to be nonaromcnt > minacc
        elseif issubset(nq,
                QueryFormula(:not, QueryFormula(:hydrogenconnected, 0)), eval_recursive=false)
            nonaromcnt += 1
        end
        # incidences
        hasarombond = false
        for inc in incidences(qmol, n)
            eq = edgeattr(qmol_, inc).query
            if issubset(eq, QueryFormula(:isaromatic, true), eval_recursive=false)
                hasarombond = true
                break
            end
            if issubset(eq, QueryFormula(:order, 2), eval_recursive=false)
                # C=O special case
                adjq = nodeattr(qmol, neighbors(qmol, n)[inc]).query
                if issubset(adjq, QueryFormula(:symbol, :O), eval_recursive=false)
                    continue
                end
            end
            if (issubset(eq, QueryFormula(:isaromatic, false), eval_recursive=false)
                    || issubset(eq, QueryFormula(:isring, false), eval_recursive=false))
                if issubset(eq, QueryFormula(:not, QueryFormula(:order, 1)), eval_recursive=false)
                    nonaromcnt += 2  # 2 is enough to be nonaromcnt > minacc
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
    aromf = QueryFormula(:isaromatic, true)
    aors = QueryFormula(:or, Set([
        aromf,
        QueryFormula(:and, Set([
            QueryFormula(:order, 1),
            QueryFormula(:isaromatic, false)
        ]))
    ]))
    aord = QueryFormula(:or, Set([
        aromf,
        QueryFormula(:and, Set([
            QueryFormula(:order, 2),
            QueryFormula(:isaromatic, false)
        ]))
    ]))
    for ring in sssr(qmol)
        ringedges = edgeset(nodesubgraph(qmol, ring))
        pcnt = 0
        for n in ring
            nq = nodeattr(qmol, n).query
            if issubset(nq, QueryFormula(:isaromatic, true), eval_recursive=false)
                pcnt += 1
                continue
            end
            rincs = collect(intersect(incidences(qmol, n), ringedges))
            uq = edgeattr(qmol, rincs[1]).query
            vq = edgeattr(qmol, rincs[2]).query
            if uq == aors && vq == aord || vq == aors && uq == aord
                if issubset(nq, QueryFormula(:or,
                            Set([QueryFormula(:symbol, a) for a in [:B, :C, :N, :P, :As]])),
                        eval_recursive=false)
                    pcnt += 1
                    continue
                end
            elseif uq == aors && vq == aors
                if issubset(nq, QueryFormula(:or,
                            Set([QueryFormula(:symbol, a) for a in [:N, :O, :P, :S, :As, :Se]])),
                        eval_recursive=false)
                    pcnt += 2
                    continue
                elseif issubset(nq, QueryFormula(:symbol, :C), eval_recursive=false)
                    outer = collect(setdiff(incidences(qmol, n), ringedges))
                    if length(outer) == 1
                        outerq = edgeattr(qmol, outer[1]).query
                        oadjq = nodeattr(qmol, neighbors(qmol, n)[outer[1]]).query
                        if issubset(outerq, QueryFormula(:order, 2), eval_recursive=false) && issubset(oadjq, QueryFormula(:symbol, :O), eval_recursive=false)
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
    # TODO: get rid of mutable Dict
    outneighbormap::Vector{Dict{Int,Int}}
    inneighbormap::Vector{Dict{Int,Int}}
    edges::Vector{Tuple{Int,Int}}
    nodeattrs::Vector{Dict}
    edgeattrs::Vector{Dict}
end
"""

"""
    dictdigraph(view::DiSubgraphView{DictDiGraph}) -> DictDiGraph

Generate a new `DictDiGraph` from a substructure view.

Graph property caches and attributes are not inherited.

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

"""
    query_relationship(;sourcefile=DEFAULT_QUERY_RELATIONS) -> DictDiGraph

Generate query relationship diagram.

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

"""
    filter_queries(qr::DictDiGraph, mol::GraphMol) -> DictDiGraph

Filter query relationship diagram by the given molecule.
The filtered diagram represents query relationship that the molecule have.

function filter_queries(qr::DictDiGraph, mol::GraphMol; filtering=true)
    matched = Set{Int}()
    for n in reversetopologicalsort(qr)
        rcd = nodeattr(qr, n)
        if filtering
            if !issubset(successors(qr, n), matched)  # query containment filter
                continue
            end
        end
        # println("key: \$(rcd["key"])")
        # println("query: \$(rcd["query"])")
        # @time begin
            matches = collect(substructmatches(mol, rcd["parsed"]))
            if !isempty(matches)
                push!(matched, n)
                rcd["matched"] = Set([sort(collect(keys(m))) for m in matches])
            end
        # end
    end
    return dictdigraph(nodesubgraph(qr, matched))
end
"""
