#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export query_containment_diagram, filter_queries


const DEFAULT_QUERIES = let
    qrfile = joinpath(dirname(@__FILE__), "../assets/const/default_queries.yaml")
    include_dependency(qrfile)
    qrfile
end



"""
    resolve_disjoint_not(tree, propmap) -> QueryTree

Resolve relationship between a `:not` query and a `:or` query and returns `QueryTree` that
the relationship is resolved.

For example, the results of `[!#7]` (not nitrogen) is included in the results of `[#8,#16]`
(oxygen or sulfur), so `[!#7]` can be converted to an equivalent query `[!#7,#8,#16]` to
generate a truthtable with common properties.
"""
function resolve_disjoint_not(tree, propmap)
    tree isa QueryOperator || return tree
    if tree.key === :not
        cld = tree.value[1]
        haskey(propmap, cld.key) || return tree
        vals = []
        for p in propmap[cld.key]
            p.value == cld.value && continue
            push!(vals, p)
        end
        return isempty(vals) ? tree : QueryOperator(:or, [tree, vals...])
    else
        return QueryOperator(tree.key,
            [resolve_disjoint_not(v, propmap) for v in tree.value])
    end
end


"""
    resolve_recursive(tree, propmap) -> QueryTree

Resolve relationship between `:recursive` queries and returns `QueryTree` that the
relationship is resolved.

For example, the results of `[\$(C)]` is included in the results of `[\$(CO)]`,
so `[\$(C)]` can be converted to an equivalent query ``[\$(C),\$(CO)]` to
generate a truthtable with common properties.
"""
function resolve_recursive(tree, propmap, rec_cache)
    tree isa QueryAny && return tree
    if tree.key === :recursive
        vals = Union{QueryAny,QueryLiteral,QueryOperator}[]
        !isnothing(rec_cache) && !haskey(rec_cache, tree.value) && (
            rec_cache[tree.value] = Dict{String,Bool}())
        tmol = smartstomol(tree.value)
        push!(vals, get_prop(tmol, 1, :tree))  # $(CN) == [$(CN);C]
        if haskey(propmap, :recursive)
            for p in propmap[:recursive]
                p.value == tree.value && continue
                if !isnothing(rec_cache) && haskey(rec_cache[tree.value], p.value)
                    match = rec_cache[tree.value][p.value]  # use cache
                else
                    match = has_substruct_match(tmol, smartstomol(p.value);
                        mandatory=Dict(1 => 1))
                    isnothing(rec_cache) || (rec_cache[tree.value][p.value] = match)  # set cache
                end
                match && push!(vals, p)
            end
        end
        return QueryOperator(:and, [tree, vals...])
    elseif tree.key === :not
        return QueryOperator(tree.key,
            [resolve_recursive(tree.value[1], propmap, rec_cache)])
    elseif tree.key in (:and, :or)
        return QueryOperator(tree.key,
            [resolve_recursive(v, propmap, rec_cache) for v in tree.value])
    else
        return tree
    end
end


"""
    generate_truthtable(q, r; recursive_cache=nothing) -> QueryTruthTable

Generate truthtable to compare queries.
"""
function generate_truthtable(q, r; recursive_cache=nothing, kwargs...)
    qt = optimize_query(q.tree)
    rt = optimize_query(r.tree)
    # convert not query (e.g. !C -> [!#6,!A]) smarts/logicaloperator.jl
    qpropmap = querypropmap(qt)
    rpropmap = querypropmap(rt)
    # resolve disjoint not (e.g. [#7,#8] => [!#6] ---> [#7,#8] => [!#6,#7,#8])
    qt = resolve_disjoint_not(qt, rpropmap)
    rt = resolve_disjoint_not(rt, qpropmap)
    # resolve recursive (e.g. $(CN) => $(C[NH2]) ---> $(CN) => [$(C[NH2]);$(CN)]
    # this updates props, so need querypropmap recalculation
    qt = resolve_recursive(qt, rpropmap, recursive_cache)
    rt = resolve_recursive(rt, qpropmap, recursive_cache)
    # reconstruct functions
    props = sort(union(
        QueryLiteral[], values(querypropmap(qt))..., values(querypropmap(rt))...))
    qfunc = generate_queryfunc(qt, props)
    rfunc = generate_queryfunc(rt, props)
    return (QueryTruthTable(qfunc, props), QueryTruthTable(rfunc, props))
end


function querymatch(q::QueryTruthTable, r::QueryTruthTable, exactmatch; maxsize=14, kwargs...)
    # truth table vector match
    # TODO: naive implementation costs worst O(2^n)
    q.props == r.props || error("query property mismatch")
    nlit = length(q.props)
    nlit > maxsize && (@info "MolecularGraph.querymatch: maxsize exceeded"; return false)
    for i in 1:(2^nlit)
        arr = isone.(digits(i-1, base=2, pad=nlit)[1:nlit])
        qout = q.func(arr)
        rout = r.func(arr)
        qout && !rout && return false
        exactmatch && !qout && rout && return false
    end
    return true
end



function querymatch(q::QueryTree, r::QueryTree, exactmatch; kwargs...)
    qtbl, rtbl = generate_truthtable(q, r; kwargs...)
    return querymatch(qtbl, rtbl, exactmatch; kwargs...)
end

"""
    Base.:(==)(q::QueryTruthTable, r::QueryTruthTable; kwargs...) 

Returns whether the two queries are equivalent
"""
Base.:(==)(q::QueryTruthTable, r::QueryTruthTable; kwargs...) = querymatch(q, r, true; kwargs...)
Base.:(==)(q::QueryTree, r::QueryTree; kwargs...) = querymatch(q, r, true; kwargs...)

"""
    Base.issubset(q::QueryTruthTable, r::QueryTruthTable; kwargs...) 

Returns whether all the results of query `q` is included in the results of query `r`
"""
Base.issubset(q::QueryTruthTable, r::QueryTruthTable; kwargs...) = querymatch(q, r, false; kwargs...)
Base.issubset(q::QueryTree, r::QueryTree; kwargs...) = querymatch(q, r, false; kwargs...)



"""
    query_containment_diagram(;sourcefile=DEFAULT_QUERY_DEFAULT_QUERIES
        ) -> DictDiGraph, vprops, eprops

Generate query containment diagram.
"""
function query_containment_diagram(;sourcefile=DEFAULT_QUERIES)
    # TODO: should be MetaGraph
    g = SimpleDiGraph{Int}()
    vs = Dict{Int,Dict}()
    es = Dict{Edge{Int},Symbol}()
    nrevmap = Dict{String,Int}()
    for (i, rcd) in enumerate(YAML.load(open(sourcefile)))
        rcd["parsed"] = smartstomol(rcd["query"])
        nrevmap[rcd["key"]] = i
        add_vertex!(g)
        vs[i] = rcd
    end
    for i in vertices(g)
        rcd = vs[i]
        if haskey(rcd, "isa")
            for k in rcd["isa"]
                e = Edge{Int}(i => nrevmap[k])
                add_edge!(g, e)
                es[e] = :isa
            end
        end
        if haskey(rcd, "has")
            for k in rcd["has"]
                e = Edge{Int}(i => nrevmap[k])
                add_edge!(g, e)
                es[e] = :has
            end
        end
    end
    return g, vs, es
end


"""
    find_queries(mol::MolGraph, query_diagram; subsets=[], filtering=true
        ) -> DictDiGraph, vprops, eprops

Find query relationship diagram by the given molecule.
The filtered diagram represents query relationship that the molecule have.
"""
function find_queries(mol::MolGraph, query_diagram; subsets=[], filtering=true)
    qr, vs, es = query_diagram
    matched = Set{Int}()
    vs_ = deepcopy(vs)
    for n in topological_sort(reverse(qr))
        rcd = vs_[n]
        if filtering
            if !issubset(outneighbors(qr, n), matched)  # query containment filter
                continue
            end
        end
        # println("key: \$(rcd["key"])")
        # println("query: \$(rcd["query"])")
        # @time begin
        matches = collect(substruct_matches(mol, rcd["parsed"]))
        if !isempty(matches)
            push!(matched, n)
            rcd["matched"] = Set([sort(collect(keys(m))) for m in matches])
        end
        # end
    end
    subg, nmap = induced_subgraph(qr, collect(matched))
    revmap = Dict(v => i for (i, v) in enumerate(nmap))
    es_ = Dict{Edge{Int},Symbol}()
    for (k, v) in es
        src(k) in nmap && dst(k) in nmap && (es_[Edge{Int}(revmap[src(k)], revmap[dst(k)])] = v)
    end
    return subg, vs_[nmap], es_
end
