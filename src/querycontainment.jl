#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#


# query_relationship, filter_queries

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


function resolve_recursive(tree, propmap, rec_cache=nothing)
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


function querymatch(q::QueryTruthTable, r::QueryTruthTable, exactmatch; maxsize=2^20, kwargs...)
    # truth table vector match
    # TODO: naive implementation costs worst O(2^n)
    q.props == r.props || throw(ErrorException("query property mismatch"))
    nlit = length(q.props)
    for i in 1:(2^nlit)
        i > maxsize && (@info "MolecularGraph.querymatch: maxsize reached"; return false)
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

Base.:(==)(q::QueryTruthTable, r::QueryTruthTable; kwargs...) = querymatch(q, r, true; kwargs...)
Base.:(==)(q::QueryTree, r::QueryTree; kwargs...) = querymatch(q, r, true; kwargs...)
Base.issubset(q::QueryTruthTable, r::QueryTruthTable; kwargs...) = querymatch(q, r, false; kwargs...)
Base.issubset(q::QueryTree, r::QueryTree; kwargs...) = querymatch(q, r, false; kwargs...)



const DEFAULT_QUERY_RELATIONS = let
    qrfile = joinpath(dirname(@__FILE__), "../../assets/const/default_query_relations.yaml")
    include_dependency(qrfile)
    qrfile
end


"""
    query_relationship(;sourcefile=DEFAULT_QUERY_RELATIONS) -> DictDiGraph

Generate query relationship diagram.
"""
function query_relationship(;sourcefile=DEFAULT_QUERY_RELATIONS)
    # TODO: attributed graph (worth adding MetaGraph dependency?)
    g = SimpleDiGraph{Int}()
    vs = Dict{String,Any}[]
    es = Dict{Edge{Int},Symbol}()
    nrevmap = Dict{String,Int}()
    for (i, rcd) in enumerate(YAML.load(open(sourcefile)))
        rcd["parsed"] = smartstomol(rcd["query"])
        nrevmap[rcd["key"]] = i
        add_node!(g)
        push!(vs, rcd)
    end
    for i in vertices(graph)
        rcd = vs[i]
        if haskey(rcd, "isa")
            for k in rcd["isa"]
                e = Edge{Int}(i => nrevmap[k])
                addedge!(g, e)
                es[e] = :isa
            end
        end
        if haskey(rcd, "has")
            for k in rcd["has"]
                e = Edge{Int}(i => nrevmap[k])
                addedge!(g, e)
                es[e] = :has
            end
        end
    end
    return g, vs, es
end


"""
    filter_queries(qr::SimpleDiGraph, vs, es, mol::MolGraph; filtering=true) -> DictDiGraph

Filter query relationship diagram by the given molecule.
The filtered diagram represents query relationship that the molecule have.
"""
function filter_queries(qr::SimpleDiGraph, vs, es, mol::MolGraph; filtering=true)
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
    subg, nmap = induced_subgraph(qr, matched)
    revmap = Dict(v => i for (i, v) in enumerate(nmap))
    es_ = Dict{Edge{Int},Symbol}()
    for (k, v) in es
        src(k) in nmap && dst(k) in nmap && (es_[Edge{Int}(revmap[src(k)], revmap[dst(k)])] = v)
    end
    return subg, vs_[nmap], es_
end
