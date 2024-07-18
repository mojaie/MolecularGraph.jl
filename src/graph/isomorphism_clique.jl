#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    ConstraintArrayMCIS, ConstraintArrayMCES,
    maximum_common_subgraph, maximum_common_edge_subgraph,
    mcis_constraints, mces_constraints


abstract type ConstraintArrayMCS{T} end

struct ConstraintArrayMCIS{T,D,V,E} <: ConstraintArrayMCS{T}
    nv::Int
    pairs::Vector{Tuple{T,T}}
    distances::Vector{D}
    vattrs1::Vector{V}
    vattrs2::Vector{V}
    eattrs::Vector{E}
end

ConstraintArrayMCIS{T,D,V,E}(data::Dict) where {T,D,V,E} = ConstraintArrayMCIS{T,D,V,E}(
    data["nv"], [(p...,) for p in data["pairs"]],
    data["distances"], data["vattrs1"], data["vattrs2"], data["eattrs"])

to_dict(arr::ConstraintArrayMCIS) = Dict(
    "nv" => arr.nv,
    "pairs" => arr.pairs,
    "distances" => arr.distances,
    "vattrs1" => arr.vattrs1,
    "vattrs2" => arr.vattrs2,
    "eattrs" => arr.eattrs
)


struct ConstraintArrayMCES{T,D,V,E} <: ConstraintArrayMCS{T}
    nv::Int
    pairs::Vector{Tuple{T,T}}
    distances::Vector{D}
    vattrs1::Vector{V}
    vattrs2::Vector{V}
    eattrs::Vector{E}
    revmap::Dict{T,Edge{T}}
    delta_edges::Vector{Tuple{Edge{T},Edge{T},Edge{T}}}
    y_edges::Vector{Tuple{Edge{T},Edge{T},Edge{T}}}
end

ConstraintArrayMCES{T,D,V,E}(data::Dict) where {T,D,V,E} = ConstraintArrayMCES{T,D,V,E}(
    data["nv"], [(p...,) for p in data["pairs"]],
    data["distances"], data["vattrs1"], data["vattrs2"], data["eattrs"],
    Dict(parse(T, i) => Edge{T}(e...) for (i, e) in data["revmap"]),
    [Tuple(Edge{T}(e...) for e in t) for t in data["delta_edges"]],
    [Tuple([Edge{T}(e...) for e in t]) for t in data["y_edges"]]
)

to_dict(arr::ConstraintArrayMCES) = Dict(
    "nv" => arr.nv,
    "pairs" => arr.pairs,
    "distances" => arr.distances,
    "vattrs1" => arr.vattrs1,
    "vattrs2" => arr.vattrs2,
    "eattrs" => arr.eattrs,
    "revmap" => Dict(i => [src(e), dst(e)] for (i, e) in arr.revmap),
    "delta_edges" => [[[src(e), dst(e)] for e in t] for t in arr.delta_edges],
    "y_edges" => [[[src(e), dst(e)] for e in t] for t in arr.y_edges]
)


function connection_constraints(::Type{U}, g::SimpleGraph{T}; kwargs...) where {T,U<:Integer}
    """TODO: too slow
    only for benchmarking and connected-MCS
    """
    pairs = Tuple{T,T}[]
    arr = U[]  # connected or not (1 bit)
    for i in vertices(g)
        for j in vertices(g)
            i < j || continue
            push!(pairs, (i, j))
            push!(arr, U(has_edge(g, i, j)))
        end
    end
    return pairs, arr
end
connection_constraints(g; kwargs...) = connection_constraints(Bool, g; kwargs...)


function shortest_path_constraints(::Type{U}, g::SimpleGraph{T};
        diameter=typemax(U), kwargs...) where {T,U<:Integer}
    pairs = Tuple{T,T}[]
    arr = U[]  # distance 1-255
    for i in vertices(g)
        for (j, d) in enumerate(gdistances(g, i))
            i < j || continue
            d > diameter && continue
            push!(pairs, (i, j))
            push!(arr, U(d))
        end
    end
    return pairs, arr
end
shortest_path_constraints(g; kwargs...) = shortest_path_constraints(UInt8, g; kwargs...)


function any_path_constraints(::Type{U}, g::SimpleGraph{T};
        diameter=sizeof(U)*8, kwargs...) where {T,U<:Integer}
    """TODO: too slow (especially in the case of MCES)
    can find more common vertices but maybe :tolerant is better
    """
    pairs = Tuple{T,T}[]
    arr = U[] # Int32 as a bit set of distances 1-32
    for i in vertices(g)
        for j in vertices(g)
            i < j || continue
            dists = noweight_all_distances(g, i, j)
            if 1 in dists  # has edge
                dists = one(U)
            else
                dists = sum(2^(d-1) for d in filter(x -> x <= diameter, dists))
            end
            dists == 0 && continue
            push!(pairs, (i, j))
            push!(arr, dists)
        end
    end
    return pairs, arr
end
any_path_constraints(g; kwargs...) = any_path_constraints(UInt32, g; kwargs...)


dist_constraint_func = Dict(
    :connection => connection_constraints,
    :shortest => shortest_path_constraints,
    :any => any_path_constraints
)


function mcis_constraints(g::SimpleGraph{T}, con;
        vmatchvec=x->0, ematchvec=(x, y)->0, kwargs...) where T
    pairs, distarr = dist_constraint_func[con](g; kwargs...)
    vatt1 = Int[]
    vatt2 = Int[]
    eatt = Int[]
    for (u, v) in pairs
        push!(vatt1, vmatchvec(u))
        push!(vatt2, vmatchvec(v))
        push!(eatt, ematchvec(u, v))
    end
    return ConstraintArrayMCIS{T,eltype(distarr),Int,Int}(
        nv(g), pairs, distarr, vatt1, vatt2, eatt)
end


function mces_constraints(g::SimpleGraph{T}, con;
        vmatchvec=x->0, ematchvec=(x, y)->0, kwargs...) where T
    lg, revmap, shared = line_graph(g)
    pairs, distarr = dist_constraint_func[con](lg; kwargs...)
    vatt1 = Int[]
    vatt2 = Int[]
    eatt = Int[]
    lgvmatchvec = lgvmatchvecgen(revmap, vmatchvec, ematchvec)
    lgematchvec = lgematchvecgen(lg, shared, vmatchvec)
    for (u, v) in pairs
        push!(vatt1, lgvmatchvec(u))
        push!(vatt2, lgvmatchvec(v))
        push!(eatt, lgematchvec(u, v))
    end
    vmatch = (u, v) -> vmatchvec(u) == vmatchvec(v)
    return ConstraintArrayMCES{T,eltype(distarr),Int,Int}(
        nv(lg), pairs, distarr, vatt1, vatt2, eatt, revmap,
        delta_edges(g, vmatch), y_edges(g, vmatch))
end


function mces_postprocess(g::ConstraintArrayMCES, h::ConstraintArrayMCES)
    return function (state)
        length(state.Q) > length(state.maxsofar) || return
        # modprod reverse mapping
        lnnodemap = Dict(div(i - 1, h.nv) + 1 => mod(i - 1, h.nv) + 1 for i in state.Q)  # lg(v) => lh(v)
        emap = Dict(g.revmap[m] => h.revmap[n] for (m, n) in lnnodemap)  # g(e) => h(e)
        # Δ-Y check
        delta_y_correction!(emap, g.delta_edges, g.y_edges, h.delta_edges, h.y_edges)
        if length(emap) > length(state.maxsofar)
            state.maxsofar = emap
        end
    end
end


function connfuncgen(modprod::SimpleGraph{T}, g::ConstraintArrayMCS{T}, h::ConstraintArrayMCS{T}) where T
    eltype(g.distances) === Bool || error(
        "topological constraint cannot be used for c-clique detection")
    conn = [Set{T}() for _ in vertices(modprod)]
    disconn = [Set{T}() for _ in vertices(modprod)]
    distmap = Dict(u_edge(T, u, v) => d for ((u, v), d) in zip(g.pairs, g.distances))
    for i in vertices(modprod)
        u = div(i - 1, h.nv) + 1
        for nbr in neighbors(modprod, i)
            v = div(nbr - 1, h.nv) + 1
            container = distmap[u_edge(T, u, v)] ? conn : disconn
            push!(container[i], nbr)
        end
    end
    return n -> (conn[n], disconn[n])
end


"""
    modular_product(g::ConstraintArrayMCS, h::ConstraintArrayMCS;
        tolerance=0, kwargs...) -> SimpleGraph

Compute modular product of the two MCS constrait arrays.
"""
function modular_product(g::ConstraintArrayMCS, h::ConstraintArrayMCS;
        tolerance=0, kwargs...)
    if eltype(g.distances) === UInt8
        pred = (x, y) -> begin
            xd = g.distances[x]
            yd = h.distances[y]
            (((xd == 1) && (yd == 1) && g.eattrs[x] == h.eattrs[y])
            || ((xd > 1) && (yd > 1)
                && (xd > yd ? xd - yd : yd - xd) <= tolerance)) # abs UInt may overflow
        end
    elseif eltype(g.distances) === UInt32
        pred = (x, y) -> (
            (g.distances[x] == 1 && h.distances[y] == 1 && g.eattrs[x] == h.eattrs[y])
            || g.distances[x] & h.distances[y] > 0)  # at least one distance match
    elseif eltype(g.distances) === Bool
        pred = (x, y) -> (
            (g.distances[x] && h.distances[y] && g.eattrs[x] == h.eattrs[y])
            || !g.distances[x] && !h.distances[y])
    end
    id(i, j) = (i - 1) * h.nv + j
    prod = SimpleGraph(g.nv * h.nv)
    @inbounds for i in 1:length(g.pairs)
        @simd for j in 1:length(h.pairs)
            if pred(i, j)
                g1, g2 = g.pairs[i]
                h1, h2 = h.pairs[j]
                if g.vattrs1[i] == h.vattrs1[j] && g.vattrs2[i] == h.vattrs2[j]
                    add_edge!(prod, id(g1, h1), id(g2, h2))
                end
                if g.vattrs1[i] == h.vattrs2[j] && g.vattrs2[i] == h.vattrs1[j]
                    add_edge!(prod, id(g1, h2), id(g2, h1))
                end
            end
        end
    end
    return prod
end



"""
    maximum_common_subgraph(g::ConstraintArrayMCIS, h::ConstraintArrayMCIS;
        connected=false, kwargs...) -> (Dict, Symbol)

Compute maximum common node-induced subgraph (MCIS) between the two MCS constraints.
"""
function maximum_common_subgraph(
        g::ConstraintArrayMCIS, h::ConstraintArrayMCIS; connected=false, kwargs...)
    prod = modular_product(g, h; kwargs...)
    if connected
        maxclique, status = maximum_conn_clique(prod, connfunc=connfuncgen(prod, g, h); kwargs...)
    else
        maxclique, status = maximum_clique(prod; kwargs...)
    end
    # modprod reverse mapping
    nmap = Dict(div(i - 1, h.nv) + 1 => mod(i - 1, h.nv) + 1 for i in maxclique)
    return nmap, status
end


"""
    maximum_common_subgraph(g::ConstraintArrayMCES, h::ConstraintArrayMCES;
        connected=false, kwargs...) -> (Dict, Symbol)

Compute maximum common edge-induced subgraph (MCES) between the two MCS constraints.
"""
function maximum_common_subgraph(
        g::ConstraintArrayMCES{T,D,V,E}, h::ConstraintArrayMCES{T,D,V,E};
        connected=false, kwargs...) where {T,D,V,E}
    prod = modular_product(g, h; kwargs...)
    if connected
        return maximum_conn_clique(
            Dict{Edge{T},Edge{T}}, prod, connfunc=connfuncgen(prod, g, h),
            postprocess=mces_postprocess(g, h); kwargs...)
    else
        return maximum_clique(
            Dict{Edge{T},Edge{T}}, prod, postprocess=mces_postprocess(g, h); kwargs...)
    end
end




############################
# prototype implementation
############################


function modprod_edge_filter(g, h, ematch; kwargs...)
    return function (g1, g2, h1, h2)
        has_edge(g, g1, g2) === has_edge(h, h1, h2) || return false
        !has_edge(g, g1, g2) && return true
        return ematch(u_edge(g, g1, g2), u_edge(h, h1, h2))
    end
end

function tp_constraint_filter(g, h, ematch; diameter=max(nv(g), nv(h)), tolerance=0, kwargs...)
    gdist = [gdistances(g, i) for i in vertices(g)]  # typemax(int) if not connected
    hdist = [gdistances(h, i) for i in vertices(h)]
    return function (g1, g2, h1, h2)
        has_edge(g, g1, g2) === has_edge(h, h1, h2) || return false
        if !has_edge(g, g1, g2)
            gdist[g1][g2] > diameter && return false
            hdist[h1][h2] > diameter && return false
            abs(gdist[g1][g2] - hdist[h1][h2]) > tolerance && return false
            return true
        end
        return ematch(u_edge(g, g1, g2), u_edge(h, h1, h2))
    end
end


function mcis_modular_product(g::SimpleGraph, h::SimpleGraph;
        vmatch=(gn,hn)->true, ematch=(ge,he)->true, topological=false, kwargs...)
    modprodfunc = topological ? tp_constraint_filter : modprod_edge_filter
    return modular_product(g, h,
        vmatch=vmatch, edgefilter=modprodfunc(g, h, ematch; kwargs...))
end


function maximum_common_subgraph(g::SimpleGraph, h::SimpleGraph; connected=false, kwargs...)
    prod, eattr = mcis_modular_product(g, h; kwargs...)
    if connected
        maxclique, status = maximum_conn_clique(prod, connfunc=connfuncgen(prod, eattr); kwargs...)
    else
        maxclique, status = maximum_clique(prod; kwargs...)
    end
    # modprod reverse mapping
    mapping = Dict(div(i - 1, nv(h)) + 1 => mod(i - 1, nv(h)) + 1 for i in maxclique)
    return mapping, status
end


function mces_modular_product(g::SimpleGraph, h::SimpleGraph;
        vmatch=(gn,hn)->true, ematch=(ge,he)->true,
        ggmatch=(n1,n2)->true, hhmatch=(n1,n2)->true, topological=false, kwargs...)
    # generate line graph
    lg, grev, gsh = line_graph(g)
    lh, hrev, hsh = line_graph(h)
    lvm = lgvmatch(lg, lh, grev, hrev, vmatch, ematch)
    lem = lgematch(lg, lh, gsh, hsh, vmatch)
    modprodfunc = topological ? tp_constraint_filter : modprod_edge_filter
    prod, eattr = modular_product(lg, lh,
        vmatch=lvm, edgefilter=modprodfunc(lg, lh, lem; kwargs...))
    gd, gy, hd, hy = (delta_edges(g, ggmatch), y_edges(g, ggmatch), delta_edges(h, hhmatch), y_edges(h, hhmatch))
    return prod, eattr, grev, hrev, gd, gy, hd, hy
end


function mces_postprocess(grev, hrev, gd, gy, hd, hy, neh)
    return function (state)
        length(state.Q) > length(state.maxsofar) || return
        # modprod reverse mapping
        lnnodemap = Dict(div(i - 1, neh) + 1 => mod(i - 1, neh) + 1 for i in state.Q)  # lg(v) => lh(v)
        emap = Dict(grev[m] => hrev[n] for (m, n) in lnnodemap)  # g(e) => h(e)
        # Δ-Y check
        delta_y_correction!(emap, gd, gy, hd, hy)
        if length(emap) > length(state.maxsofar)
            state.maxsofar = emap
        end
    end
end


function maximum_common_edge_subgraph(
        g::SimpleGraph{T}, h::SimpleGraph{T}; connected=false, kwargs...) where T
    prod, eattr, grev, hrev, gd, gy, hd, hy = mces_modular_product(g, h; kwargs...)
    pp = mces_postprocess(grev, hrev, gd, gy, hd, hy, ne(h))
    if connected
        return maximum_conn_clique(Dict{Edge{T},Edge{T}}, prod, connfunc=connfuncgen(prod, eattr), postprocess=pp; kwargs...)
    else
        return maximum_clique(Dict{Edge{T},Edge{T}}, prod, postprocess=pp; kwargs...)
    end
end
