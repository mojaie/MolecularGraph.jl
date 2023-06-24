#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    maximum_common_subgraph, maximum_common_edge_subgraph,
    mcis_constraints, mces_constraints


struct MCSResult{T}
    mapping::Dict{T,T}
    status::Symbol
end

Base.size(result::MCSResult) = length(result.mapping)


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


function mcs_clique_detection(prod, isconn; method=:disconn, kwargs...)
    if method === :disconn
        res = find_cliques(prod; kwargs...)
        return res.cliques, res.status
    elseif method === :connected
        res = find_conn_cliques(prod, isconn; kwargs...)
        return res.cliques, res.status
    elseif method === :approx
        return approx_maximal_cliques(prod), :done
    end
end


"""
    mcs_modular_product(g::SimpleGraph, h::SimpleGraph; kwargs...) -> SimpleGraph

Compute modular product for MCIS.
"""
function mcis_modular_product(g::SimpleGraph, h::SimpleGraph;
        vmatch=(gn,hn)->true, ematch=(ge,he)->true, topological=false, kwargs...)
    modprodfunc = topological ? tp_constraint_filter : modprod_edge_filter
    return modular_product(g, h,
        vmatch=vmatch, edgefilter=modprodfunc(g, h, ematch; kwargs...))
end


"""
    maximum_common_subgraph(g::SimpleGraph, h::SimpleGraph; kwargs...) -> MCSResult

Compute maximum common node-induced subgraph (MCIS) between graphs `g` and `h`.
"""
function maximum_common_subgraph(g::SimpleGraph, h::SimpleGraph; kwargs...)
    prod, isconn = mcis_modular_product(g, h; kwargs...)
    return _maximum_common_subgraph(prod, isconn, nv(h); kwargs...)
end

function _maximum_common_subgraph(prod, isconn, nvh; kwargs...)
    cliques, status = mcs_clique_detection(prod, isconn; kwargs...)
    maxclique = sortstablemax(cliques, by=length, init=[])
    # modprod reverse mapping
    mapping = Dict(div(i - 1, nvh) + 1 => mod(i - 1, nvh) + 1 for i in maxclique)
    return MCSResult(mapping, status)
end


"""
    maximum_common_edge_subgraph(g::SimpleGraph, h::SimpleGraph; kwargs...) -> MCSResult

Compute maximum common edge-induced subgraph (MCES) between graphs `g` and `h`.
"""
function mces_modular_product(g::SimpleGraph, h::SimpleGraph;
        vmatch=(gn,hn)->true, ematch=(ge,he)->true,
        ggmatch=(n1,n2)->true, hhmatch=(n1,n2)->true, topological=false, kwargs...)
    # generate line graph
    lg, grev, gsh = line_graph(g)
    lh, hrev, hsh = line_graph(h)
    lvm = lgvmatch(lg, lh, grev, hrev, vmatch, ematch)
    lem = lgematch(lg, lh, gsh, hsh, vmatch)
    modprodfunc = topological ? tp_constraint_filter : modprod_edge_filter
    prod, isconn = modular_product(lg, lh,
        vmatch=lvm, edgefilter=modprodfunc(lg, lh, lem; kwargs...))
    gd, gy, hd, hy = (delta_edges(g, ggmatch), y_edges(g, ggmatch), delta_edges(h, hhmatch), y_edges(h, hhmatch))
    return prod, isconn, grev, hrev, gd, gy, hd, hy
end


"""
    maximum_common_edge_subgraph(g::SimpleGraph, h::SimpleGraph; kwargs...) -> MCSResult

Compute maximum common edge-induced subgraph (MCES) between graphs `g` and `h`.
"""
function maximum_common_edge_subgraph(
        g::SimpleGraph{T}, h::SimpleGraph{T}; kwargs...) where T
    prod, isconn, grev, hrev, gd, gy, hd, hy = mces_modular_product(g, h; kwargs...)
    return _maximum_common_edge_subgraph(g, h, prod, isconn, grev, hrev, gd, gy, hd, hy, ne(h); kwargs...)
end

function _maximum_common_edge_subgraph(g::SimpleGraph{T}, h::SimpleGraph{T},
        prod, isconn, grev, hrev, gd, gy, hd, hy, neh; kwargs...) where T
    cliques, status = mcs_clique_detection(prod, isconn; kwargs...)
    mapping = Dict{Edge{T},Edge{T}}()
    for cq in cliques
        length(cq) > length(mapping) || continue
        # modprod reverse mapping
        lnnodemap = Dict(div(i - 1, neh) + 1 => mod(i - 1, neh) + 1 for i in cq)  # lg(v) => lh(v)
        emap = Dict(grev[m] => hrev[n] for (m, n) in lnnodemap)  # g(e) => h(e)
        # Δ-Y check
        delta_y_correction!(emap, gd, gy, hd, hy)
        if length(emap) > length(mapping)
            mapping = emap
        end
    end
    return MCSResult(mapping, status)
end



abstract type ConstraintArrayMCS{T} end

struct ConstraintArrayMCIS{T} <: ConstraintArrayMCS{T}
    nv::Int
    pairs::Vector{Tuple{T,T}}
    distances::Vector{BitVector}
    attributes::Vector{BitVector}
end

struct ConstraintArrayMCES{T} <: ConstraintArrayMCS{T}
    nv::Int
    pairs::Vector{Tuple{T,T}}
    distances::Vector{BitVector}
    attributes::Vector{BitVector}
    revmap::Dict{T,Edge{T}}
    delta_edges::Vector{Tuple{Edge{T},Edge{T},Edge{T}}}
    y_edges::Vector{Tuple{Edge{T},Edge{T},Edge{T}}}
end


function connection_constraints(g::SimpleGraph{T}) where T
    """TODO: too slow
    only for benchmarking and connected-MCS
    """
    pairs = Tuple{T,T}[]
    arr = BitVector[]  # connected or not (1 bit)
    for i in vertices(g)
        for j in vertices(g)
            i == j && continue
            push!(pairs, (i, j))
            push!(arr, BitVector([has_edge(g, i, j)]))
        end
    end
    return pairs, arr
end

function shortest_path_constraints(g::SimpleGraph{T};
        diameter=typemax(UInt8)+1) where T
    pairs = Tuple{T,T}[]
    arr = BitVector[]  # distance 1-256 (8 bits)
    for i in vertices(g)
        for (j, d) in enumerate(gdistances(g, i))
            i == j && continue
            d > diameter && continue
            push!(pairs, (i, j))
            push!(arr, BitVector(digits(d - 1, base=2, pad=8)))
        end
    end
    return pairs, arr
end

function tolerant_path_constraints(g::SimpleGraph{T};
        diameter=32, tolerance=1) where T
    """TODO: significantly slower than :shortest
    but can find more common vertices/edges
    """
    pairs = Tuple{T,T}[]
    arr = BitVector[] # distance bit set 1-32 (32 bits)
    for i in vertices(g)
        for (j, d) in enumerate(gdistances(g, i))
            i == j && continue
            if d == 1  # has edge
                dists = [1]
            else
                dists = [((d - tolerance):(d - 1))..., (d:(d + tolerance))...]
                dists = filter(x -> x > 1 && x <= diameter, dists)
            end
            isempty(dists) && continue
            vec = BitVector(digits(0, base=2, pad=32))
            vec[dists] .= 1
            push!(pairs, (i, j))
            push!(arr, vec)
        end
    end
    return pairs, arr
end

function any_path_constraints(g::SimpleGraph{T}; diameter=32) where T
    """TODO: too slow (especially in the case of MCES)
    can find more common vertices but maybe :tolerant is better
    """
    pairs = Tuple{T,T}[]
    arr = BitVector[] # distance bit set 1-32 (32 bits)
    for i in vertices(g)
        for j in vertices(g)
            i == j && continue
            dists = noweight_all_distances(g, i, j)
            if 1 in dists  # has edge
                dists = [1]
            else
                dists = filter(x -> x <= diameter, dists)
            end
            isempty(dists) && continue
            vec = BitVector(digits(0, base=2, pad=32))
            vec[dists] .= 1
            push!(pairs, (i, j))
            push!(arr, vec)
        end
    end
    return pairs, arr
end


dist_constraint_func = Dict(
    :connection => connection_constraints,
    :shortest => shortest_path_constraints,
    :tolerant => tolerant_path_constraints,
    :any => any_path_constraints
)


function mcis_constraints(g::SimpleGraph{T}, con;
        vmatchvec=v->BitVector(), kwargs...) where T
    # TODO: ematchvec?
    pairs, distarr = dist_constraint_func[con](g; kwargs...)
    arr = BitVector[] # atomnumber 0-255 + pi_electron 0-3 x2 (20 bits)
    for (u, v) in pairs
        push!(arr, vcat(vmatchvec(u), vmatchvec(v)))
    end
    return ConstraintArrayMCIS{T}(nv(g), pairs, distarr, arr)
end


function mces_constraints(g::SimpleGraph{T}, con;
        vmatchvec=v->BitVector(), kwargs...) where T
    # TODO: ematchvec?
    lg, revmap, shared = line_graph(g)
    pairs, distarr = dist_constraint_func[con](lg; kwargs...)
    arr = BitVector[] # size: vmatchvec x5
    lgvmatchvec = lgvmatchvecgen(revmap, vmatchvec)
    lgematchvec = lgematchvecgen(lg, shared, vmatchvec)
    for (i, j) in pairs
        vec1 = lgvmatchvec(i)
        vec2 = lgvmatchvec(j)
        svec = lgematchvec(i, j)
        push!(arr, vcat(vec1, vec2, svec))
    end
    vmatch = (u, v) -> vmatchvec(u) == vmatchvec(v)
    return ConstraintArrayMCES{T}(nv(lg), pairs, distarr, arr, revmap,
        delta_edges(g, vmatch), y_edges(g, vmatch))
end


function modular_product(g::ConstraintArrayMCS, h::ConstraintArrayMCS)
    prod = SimpleGraph(g.nv * h.nv)
    isempty(g.distances) && return prod
    if length(g.distances[1]) <= 8  # single dist
        # attribute+distance bits match
        gda = vcat.(g.distances, g.attributes)
        hda = reshape(vcat.(h.distances, h.attributes), 1, :)
        edges = gda .== hda
    else  # multi dist
        # attribute bits match
        attrmat = g.attributes .== reshape(h.attributes, 1, :)
        # AND any of the possible distance bits match
        distmat = any.(broadcast(.&, g.distances, reshape(h.distances, 1, :)))
        edges = attrmat .& distmat
    end
    id(i, j) = (i - 1) * h.nv + j
    for ci in findall(edges)
        g1, g2 = g.pairs[ci[1]]
        h1, h2 = h.pairs[ci[2]]
        add_edge!(prod, id(g1, h1), id(g2, h2))
    end
    return prod
end


function maximum_common_subgraph(
        g::ConstraintArrayMCIS, h::ConstraintArrayMCIS; kwargs...)
    prod = modular_product(g, h)
    state = find_cliques(prod; kwargs...)
    maxclique = sortstablemax(state.cliques, by=length, init=[])
    # modprod reverse mapping
    mapping = Dict(div(i - 1, h.nv) + 1 => mod(i - 1, h.nv) + 1 for i in maxclique)
    return MCSResult(mapping, state.status)
end


function maximum_common_subgraph(
        g::ConstraintArrayMCES{T}, h::ConstraintArrayMCES{T}; kwargs...) where T
    prod = modular_product(g, h)
    state = find_cliques(prod; kwargs...)
    mapping = Dict{Edge{T},Edge{T}}()
    for cq in state.cliques
        length(cq) > length(mapping) || continue
        # modprod reverse mapping
        lnnodemap = Dict(div(i - 1, h.nv) + 1 => mod(i - 1, h.nv) + 1 for i in cq)  # lg(v) => lh(v)
        emap = Dict(g.revmap[m] => h.revmap[n] for (m, n) in lnnodemap)  # g(e) => h(e)
        # Δ-Y check
        delta_y_correction!(emap, g.delta_edges, g.y_edges, h.delta_edges, h.y_edges)
        if length(emap) > length(mapping)
            mapping = emap
        end
    end
    return MCSResult(mapping, state.status)
end
