#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

abstract type ConstraintArrayMCS{T} end

struct ConstraintArrayMCIS{T,D,V,E} <: ConstraintArrayMCS{T}  # Serializable
    nv::Int
    pairs::Vector{Tuple{T,T}}
    distances::Vector{D}
    vattrs1::Vector{V}
    vattrs2::Vector{V}
    eattrs::Vector{E}
end


struct ConstraintArrayMCES{T,D,V,E} <: ConstraintArrayMCS{T}  # Serializable
    nv::Int
    pairs::Vector{Tuple{T,T}}
    distances::Vector{D}
    vattrs1::Vector{V}
    vattrs2::Vector{V}
    eattrs::Vector{E}
    revmap::Dict{VertexKey{T},Edge{T}}
    delta_edges::Vector{Tuple{Edge{T},Edge{T},Edge{T}}}
    y_edges::Vector{Tuple{Edge{T},Edge{T},Edge{T}}}
end


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


function mcis_constraints(pairs, distarr::Vector{T}, g::SimpleGraph{U};
        vmatchvec=x->0, ematchvec=(x, y)->0, kwargs...) where {T,U}
    vatt1 = Int[]
    vatt2 = Int[]
    eatt = Int[]
    for (u, v) in pairs
        push!(vatt1, vmatchvec(u))
        push!(vatt2, vmatchvec(v))
        push!(eatt, ematchvec(u, v))
    end
    return ConstraintArrayMCIS{U,T,Int,Int}(
        nv(g), pairs, distarr, vatt1, vatt2, eatt)
end

function mcis_constraints(::Val{:connection}, g::SimpleGraph; kwargs...)
    pairs, distarr = connection_constraints(g; kwargs...)
    return mcis_constraints(pairs, distarr, g; kwargs...)
end

function mcis_constraints(::Val{:shortest}, g::SimpleGraph; kwargs...)
    pairs, distarr = shortest_path_constraints(g; kwargs...)
    return mcis_constraints(pairs, distarr, g; kwargs...)
end

function mcis_constraints(::Val{:any}, g::SimpleGraph; kwargs...)
    pairs, distarr = any_path_constraints(g; kwargs...)
    return mcis_constraints(pairs, distarr, g; kwargs...)
end


function mces_constraints(pairs, distarr::Vector{T}, lg, revmap, shared, g::SimpleGraph{U};
        vmatchvec=x->0, ematchvec=(x, y)->0, kwargs...) where {T,U}
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
    return ConstraintArrayMCES{U,T,Int,Int}(
        nv(lg), pairs, distarr, vatt1, vatt2, eatt, revmap,
        delta_edges(g, vmatch), y_edges(g, vmatch))
end

function mces_constraints(::Val{:connection}, g::SimpleGraph; kwargs...)
    lg, revmap, shared = line_graph(g)
    pairs, distarr = connection_constraints(lg; kwargs...)
    return mces_constraints(pairs, distarr, lg, revmap, shared, g; kwargs...)
end

function mces_constraints(::Val{:shortest}, g::SimpleGraph; kwargs...)
    lg, revmap, shared = line_graph(g)
    pairs, distarr = shortest_path_constraints(lg; kwargs...)
    return mces_constraints(pairs, distarr, lg, revmap, shared, g; kwargs...)
end

function mces_constraints(::Val{:any}, g::SimpleGraph; kwargs...)
    lg, revmap, shared = line_graph(g)
    pairs, distarr = any_path_constraints(lg; kwargs...)
    return mces_constraints(pairs, distarr, lg, revmap, shared, g; kwargs...)
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
