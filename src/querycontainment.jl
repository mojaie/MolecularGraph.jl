#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@kwdef mutable struct QueryTerm
    key::String
    name::String
    aliases::Vector{String} = String[]
    isa::Vector{String} = String[]
    has::Vector{String} = String[]
    source::String
    query::String
end

Base.getindex(iso::QueryTerm, k::String) = getproperty(iso, Symbol(k))

const DEFAULT_QUERIES = let
    qrfile = joinpath(dirname(@__FILE__), "../assets/const/default_queries.yaml")
    include_dependency(qrfile)
    terms = QueryTerm[]
    for rcd in YAML.load(open(qrfile))
        t = QueryTerm(; NamedTuple((Symbol(k), v) for (k, v) in rcd)...)
        push!(terms, t)
    end
    terms
end


"""
    resolve_disjoint_not!(
        tree::QueryTree, propmap::Dict{Symbol,Vector{QueryNode}}) -> Nothing

Resolve relationship between a `:not` query and a `:or` query and returns `QueryTree` that
the relationship is resolved.

For example, the results of `[!#7]` (not nitrogen) is included in the results of `[#8,#16]`
(oxygen or sulfur), so `[!#7]` can be converted to an equivalent query `[!#7,#8,#16]` to
generate a truthtable with common properties.
"""
function resolve_disjoint_not!(qtree::QueryTree, propmap::Dict{Symbol,Vector{QueryNode}})
    for node in topological_sort(reverse(qtree.graph))
        nprop = qtree.vprops[node]
        nprop.operator === :not || continue
        succ = only(outneighbors(qtree.graph, node))
        sp = qtree.vprops[succ]
        haskey(propmap, sp.key) || continue
        # bypass pred->qor->qnot
        set_qnode!(qtree, node, qor())
        rem_qedge!(qtree, node, succ)
        qn = add_qnode!(qtree, node, qnot())
        add_qedge!(qtree, qn, succ)
        for pm in propmap[sp.key]
            pm == sp && continue
            add_qnode!(qtree, node, pm)
        end
    end
    return
end


"""
    resolve_recursive(
        tree::QueryTree, propmap::Dict{Int,Vector{QueryNode}}
        ;rec_cache=Dict{String,Dict{String,Bool}}()
    ) -> Nothing

Resolve relationship between `:recursive` queries and returns `QueryTree` that the
relationship is resolved.

For example, the results of `[\$(C)]` is included in the results of `[\$(CO)]`,
so `[\$(C)]` can be converted to an equivalent query ``[\$(C);\$(CO)]` to
generate a truthtable with common properties.
"""
resolve_recursive!(
    qtree::QueryTree, propmap::Dict{Symbol,Vector{QueryNode}}
    ; rec_cache=Dict{String,Dict{String,Bool}}()
) = resolve_recursive!(qtree, propmap, rec_cache)

# `rec_cache` is positional because a keyword method inside the recursion
# (querymatch -> resolve_recursive! -> has_substruct_match -> ... -> querymatch)
# is not compiled by `juliac --trim`.
function resolve_recursive!(
        qtree::QueryTree, propmap::Dict{Symbol,Vector{QueryNode}},
        rec_cache::Dict{String,Dict{String,Bool}})
    for node in topological_sort(reverse(qtree.graph))
        nprop = qtree.vprops[node]
        nprop.key === :recursive || continue
        # insert :and, reconnect pred->:and->node
        pred = inneighbors(qtree.graph, node)
        qa = isempty(pred) ? add_qnode!(qtree, qand()) : add_qnode!(qtree, only(pred), qand())
        isempty(pred) || rem_qedge!(qtree, only(pred), node)
        add_qedge!(qtree, qa, node)
        # [$([nH]C=O)], [nH] -> [$([nH]C);[nH]]
        gnv = nv(qtree.graph)
        tmol = smartstomol(nprop.value)
        head_ = tmol[1]
        headroot = root(head_) + gnv
        disjoint_union!(qtree.graph, head_.graph)
        for (k, v) in head_.vprops
            qtree.vprops[k + gnv] = v
        end
        add_qedge!(qtree, qa, headroot)
        # [$([nH]C=O);[nH]], [$([nH]C)] -> [$([nH]C=O);[nH];$([nH]C)]
        haskey(propmap, :recursive) || continue
        if !haskey(rec_cache, nprop.value)
            rec_cache[nprop.value] = Dict{String,Bool}()
        end
        for pm in propmap[:recursive]
            pm.value == nprop.value && continue
            if haskey(rec_cache[nprop.value], pm.value)
                match = rec_cache[nprop.value][pm.value]  # use cache
            else
                match = has_substruct_match(tmol, smartstomol(pm.value);
                    mandatory=Dict(1 => 1))
                rec_cache[nprop.value][pm.value] = match  # set cache
            end
            match && add_qnode!(qtree, qa, pm)
        end
    end
    return
end


function to_cnf(
        g::QueryTree, varids::Dict{Int,Int}, nextid::Int;
        inverse::Bool=false, assert_root::Bool=true)
    # Plaisted–Greenbaum(PG) transformation
    clauses = Vector{Int}[]
    nodeids = Dict{Int,Int}()

    # Post-order traversal
    # state: 0 = unvisited, 1 = visiting, 2 = done
    state = zeros(UInt8, nv(g.graph))
    order = Int[]
    r = root(g)
    stack = Tuple{Int,Bool}[(r, false)]
    while !isempty(stack)
        v, expanded = pop!(stack)
        if expanded
            state[v] = 2
            push!(order, v)
            continue
        end
        state[v] == 2 && continue
        state[v] = 1
        push!(stack, (v, true))
        for u in reverse(outneighbors(g.graph, v))
            push!(stack, (u, false))
        end
    end

    # Assign literals to nodes
    for v in order
        children = outneighbors(g.graph, v)
        if g.vprops[v].operator == :eq
            if g.vprops[v].key in (:stereo, :direction)
                # ignore stereo
                nodeids[v] = nextid
                push!(clauses, [nextid])
                nextid += 1
            else
                nodeids[v] = varids[v]
            end
        elseif g.vprops[v].operator == :any
            nodeids[v] = nextid
            push!(clauses, [nextid])
            nextid += 1
        elseif g.vprops[v].operator == :not
            # Reuse a negated child literal.
            nodeids[v] = -nodeids[only(children)]
        elseif g.vprops[v].operator in (:and, :or)
            nodeids[v] = nextid
            nextid += 1
        end
    end

    # Polarity-aware PG encoding
    sign = inverse ? -1 : 1
    # (vertex, polarity)
    stack = Tuple{Int,Int}[(r, sign)]
    visited = Set{Tuple{Int,Int}}()
    while !isempty(stack)
        v, pol = pop!(stack)
        (v, pol) in visited && continue
        push!(visited, (v, pol))
        g.vprops[v].operator in (:eq, :any) && continue
        children = outneighbors(g.graph, v)
        if g.vprops[v].operator == :not
            # Negation reverses polarity
            push!(stack, (only(children), -pol))
            continue
        end
        x = nodeids[v]
        lits = [nodeids[u] for u in children]
        if g.vprops[v].operator == :and
            if pol == 1
                # x => AND(lits)
                for c in lits
                    push!(clauses, [-x, c])
                end
            else
                # AND(lits) => x
                push!(clauses, [x; -lits])
            end
        elseif g.vprops[v].operator == :or
            if pol == 1
                # x => OR(lits)
                push!(clauses, [-x; lits])
            else
                # OR(lits) => x
                for c in lits
                    push!(clauses, [x, -c])
                end
            end
        end
        # Propagate polarity to children.
        for u in children
            push!(stack, (u, pol))
        end
    end

    # Step 4: Assert the root
    rootlit = nodeids[r]
    if assert_root
        push!(clauses, [sign * rootlit])
    end
    return clauses, nextid
end


function querymatch(
        q::T, r::U, exactmatch::Bool
        ; recursive_cache=Dict{String,Dict{String,Bool}}()) where {T<:QueryTree,U<:QueryTree}
    q = copy(q)
    r = copy(r)
    # convert not query (e.g. !C -> [!#6,!A]) smarts/logicaloperator.jl
    qpropmap = querypropmap(q)
    rpropmap = querypropmap(r)
    # resolve disjoint not (e.g. [#7,#8] => [!#6] ---> [#7,#8] => [!#6,#7,#8])
    resolve_disjoint_not!(q, rpropmap)
    resolve_disjoint_not!(r, qpropmap)
    # resolve recursive (e.g. $(CN) => $(C[NH2]) ---> $(CN) => [$(C[NH2]);$(CN)]
    # this updates props, so need querypropmap recalculation
    resolve_recursive!(q, rpropmap, recursive_cache)
    resolve_recursive!(r, qpropmap, recursive_cache)
    varid_map = Dict{QueryNode,Int}()  # Query -> id
    varids_q = Dict{Int,Int}()  # vertex -> id
    nextid = 1
    for i in vertices(q.graph)
        node = q.vprops[i]
        node.operator == :eq || continue
        if haskey(varid_map, node)
            varids_q[i] = varid_map[node]
        else
            varids_q[i] = nextid
            varid_map[node] = varids_q[i]
            nextid += 1
        end
    end
    varids_r = Dict{Int,Int}()
    for i in vertices(r.graph)
        node = r.vprops[i]
        node.operator == :eq || continue
        if haskey(varid_map, node)
            varids_r[i] = varid_map[node]
        else
            varids_r[i] = nextid
            varid_map[node] = varids_r[i]
            nextid += 1
        end
    end
    qcnf, qnext = to_cnf(q, varids_q, nextid)
    rcnf, _ = to_cnf(r, varids_r, qnext; inverse=true)
    qrcnf = [qcnf; rcnf]
    qrsat = PicoSAT.solve(qrcnf) == :unsatisfiable
    if exactmatch
        rcnf, rnext = to_cnf(r, varids_r, nextid)
        qcnf, _ = to_cnf(q, varids_q, rnext; inverse=true)
        rqcnf = [rcnf; qcnf]
        rqsat = PicoSAT.solve(rqcnf) == :unsatisfiable
        return qrsat && rqsat
    else
        return qrsat
    end
end


"""
    Base.:(==)(q::QueryTruthTable, r::QueryTruthTable; kwargs...) 

Returns whether the two queries are equivalent
"""
# Base.:(==)(q::QueryTree, r::QueryTree; kwargs...) = querymatch(q, r, true; kwargs...)


"""
    Base.issubset(q::QueryTruthTable, r::QueryTruthTable; kwargs...) 

Returns whether all the results of query `q` is included in the results of query `r`
"""
Base.issubset(q::QueryTree, r::QueryTree; kwargs...) = querymatch(q, r, false; kwargs...)


"""
    query_containment_diagram(;sourcefile=DEFAULT_QUERIES
        ) -> DictDiGraph, vprops, eprops

Generate query containment diagram.
"""
function query_containment_diagram(;sources=[], termset=DEFAULT_QUERIES)
    # filter sources
    filtered = Dict{String,QueryTerm}()
    for rcd in termset
        (isempty(sources) || (rcd.source in sources)) || continue
        filtered[rcd.key] = rcd
    end
    # merge duplicate entries
    dupes = Set{String}()
    merged = Dict{String,Dict}()
    for rcd in values(filtered)
        rcd.key in dupes && continue
        newrcd = Dict(
            "key" => rcd.key,
            "parsed" => smartstomol(rcd.query),
            "sources" => [rcd.source],
            "info" => [Dict(
                "name" => rcd.name,
                "query" => rcd.query,
                "source" => rcd.source
            )]
        )
        hasproperty(rcd, :isa) && (newrcd["isa"] = rcd.isa)
        hasproperty(rcd, :has) && (newrcd["has"] = rcd.has)
        if hasproperty(rcd, :aliases)
            for aliase in rcd.aliases
                if aliase in keys(filtered)
                    arcd = filtered[aliase]
                    push!(newrcd["info"], Dict(
                        "name" => arcd.name,
                        "query" => arcd.query,
                        "source" => arcd.source
                    ))
                end
                push!(dupes, aliase)
            end
        end
        merged[rcd.key] = newrcd
    end
    # generate graph
    # TODO: should be MetaGraph
    g = SimpleDiGraph{Int}()
    vprops = Dict{Int,Dict}()
    eprops = Dict{Edge{Int},Symbol}()
    nrevmap = Dict{String,Int}()
    for (i, (k, rcd)) in enumerate(merged)
        nrevmap[k] = i
        add_vertex!(g)
        vprops[i] = rcd
    end
    for i in vertices(g)
        if haskey(vprops[i], "isa")
            for k in vprops[i]["isa"]
                haskey(nrevmap, k) || continue  # in filtered source
                e = Edge{Int}(i => nrevmap[k])
                add_edge!(g, e)
                eprops[e] = :isa
            end
        end
        if haskey(vprops[i], "has")
            for k in vprops[i]["has"]
                haskey(nrevmap, k) || continue  # in filtered source
                e = Edge{Int}(i => nrevmap[k])
                add_edge!(g, e)
                eprops[e] = :has
            end
        end
    end
    # transitive reduction
    new_g = transitivereduction(g)
    new_eprops = Dict{Edge{Int},Symbol}()
    for e in edges(new_g)
        new_eprops[e] = eprops[e]
    end
    return new_g, vprops, new_eprops
end


"""
    find_queries(mol::MolGraph, query_diagram; subsets=[], filtering=true
        ) -> DictDiGraph, vprops, eprops

Find query relationship diagram by the given molecule.
The filtered diagram represents query relationship that the molecule have.
"""
function find_queries(mol::SimpleMolGraph, query_diagram; sources=[], filtering=true)
    qr, vs, es = query_diagram
    matched = Set{Int}()
    vs_ = copy(vs)
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
            rcd["matched"] = [collect(s) for s in Set(keys(m) for m in matches)]
        end
        # end
    end
    filtered = Int[]
    for m in matched
        (isempty(sources) || !isempty(intersect(vs_[m]["sources"], sources))) || continue
        push!(filtered, m)
    end
    subg, nmap = induced_subgraph(qr, filtered)
    revmap = Dict(v => i for (i, v) in enumerate(nmap))
    es_ = Dict{Edge{Int},Symbol}()
    for (k, v) in es
        src(k) in nmap && dst(k) in nmap && (es_[Edge{Int}(revmap[src(k)], revmap[dst(k)])] = v)
    end
    new_vs = Dict{Int,Dict}()
    for (i, v) in enumerate(nmap)
        new_vs[i] = vs_[nmap[i]]
    end
    return subg, new_vs, es_
end
