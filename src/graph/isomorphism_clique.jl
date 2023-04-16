#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    MCSResult, maximum_common_subgraph, maximum_common_edge_subgraph


struct MCSResult{T}
    mapping::Dict{T,T}
    status::Symbol
end

Base.size(result::MCSResult) = length(result.mapping)


function triangle_nodes(g::SimpleGraph{T}) where T
    ts = Set{Tuple{T,T,T}}()
    for n in vertices(g)
        nbrs = neighbors(g, n)
        length(nbrs) < 2 && continue
        for (e1, e2) in combinations(length(nbrs))
            u, v = (nbrs[e1], nbrs[e2])
            if has_edge(g, u, v)
                push!(ts, tuple(sort([n, u, v])...))
            end
        end
    end
    return ts
end


function delta_y_correction!(mapping, g, h)  # edge mapping g(e) => h(e)
    """
    any node attr mismatch would break isomorphism of triangles (cannot be matched)
    any edge attr mismatch would not affect isomorphism of other two edges in the triangle
    -> any edge can be safely removed to correct delta-Y transformation
    """
    revmap = Dict(v => k for (k, v) in mapping)
    for t in triangle_nodes(g)
        te = Edge{eltype(g)}.([(t[1], t[2]), (t[2], t[3]), (t[1], t[3])])
        issubset(te, keys(mapping)) || continue
        hs = Set{eltype(g)}()
        for e in te
            he = mapping[e]
            push!(hs, src(he), dst(he))
        end
        if length(hs) != 3
            delete!(mapping, te[1])  # remove arbitrary one edge in the triangle
        end
    end
    for t in triangle_nodes(h)
        te = Edge{eltype(h)}.([(t[1], t[2]), (t[2], t[3]), (t[1], t[3])])
        issubset(te, values(mapping)) || continue
        gs = Set{eltype(h)}()
        for e in te
            ge = revmap[e]
            push!(gs, src(ge), dst(ge))
        end
        if length(gs) != 3
            delete!(mapping, revmap[te[1]])  # remove arbitrary one edge in the triangle
        end
    end
end


function lgvmatch(g, h, grev, hrev, vmatch, ematch)
    return function (lgn, lhn)
        ge = grev[lgn]
        he = hrev[lhn]
        ematch(ge, he) || return false
        return (
            (vmatch(src(ge), src(he)) && vmatch(dst(ge), dst(he)))
            || (vmatch(src(ge), dst(he)) && vmatch(dst(ge), src(he)))
        )
    end
end


function lgematch(g, h, gsh, hsh, vmatch)
    return (lge, lhe) -> vmatch(gsh[lge], hsh[lhe])
end


function modprod_edge_filter(g, h, ematch; kwargs...)
    return function (g1, g2, h1, h2)
        has_edge(g, g1, g2) === has_edge(h, h1, h2) || return false
        !has_edge(g, g1, g2) && return true
        return ematch(undirectededge(g, g1, g2), undirectededge(h, h1, h2))
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
        return ematch(undirectededge(g, g1, g2), undirectededge(h, h1, h2))
    end
end


"""
    maximum_common_subgraph(g::SimpleGraph, h::SimpleGraph; kwargs...) -> MCSResult

Compute maximum common node-induced subgraph (MCIS) between graphs `g` and `h`.
"""
function maximum_common_subgraph(g::SimpleGraph{T}, h::SimpleGraph{T};
        vmatch=(gn,hn)->true, ematch=(ge,he)->true,
        topological=false, connected=false, kwargs...) where T
    (nv(g) == 0 || nv(h) == 0) && return MCSResult(Dict{T,T}(), :done)
    # Generate modular product
    modprodfunc = topological ? tp_constraint_filter : modprod_edge_filter
    prod, isconn = modular_product(g, h,
        vmatch=vmatch, edgefilter=modprodfunc(g, h, ematch; kwargs...))
    # Clique detection
    cqstate = (connected ?
        find_conn_cliques(prod, isconn; kwargs...) : find_cliques(prod; kwargs...))
    maxclique = sortstablemax(cqstate.cliques, by=length, init=[])
    # modprod reverse mapping
    mapping = Dict(div(i - 1, nv(h)) + 1 => mod(i - 1, nv(h)) + 1 for i in maxclique)
    return MCSResult(mapping, cqstate.status)
end


"""
    maximum_common_edge_subgraph(g::SimpleGraph, h::SimpleGraph; kwargs...) -> MCSResult

Compute maximum common edge-induced subgraph (MCES) between graphs `g` and `h`.
"""
function maximum_common_edge_subgraph(g::SimpleGraph{T}, h::SimpleGraph{T};
        vmatch=(gn,hn)->true, ematch=(ge,he)->true,
        topological=false, connected=false, kwargs...) where T
    (ne(g) == 0 || ne(h) == 0) && return MCSResult(Dict{Edge{T},Edge{T}}(), :done)
    # generate line graph
    lg, grev, gsh = line_graph(g)
    lh, hrev, hsh = line_graph(h)
    lvm = lgvmatch(lg, lh, grev, hrev, vmatch, ematch)
    lem = lgematch(lg, lh, gsh, hsh, vmatch)
    # Generate modular product
    modprodfilter = topological ? tp_constraint_filter : modprod_edge_filter
    prod, isconn = modular_product(lg, lh,
        vmatch=lvm, edgefilter=modprodfilter(lg, lh, lem; kwargs...))
    # Clique detection
    cqstate = (connected ?
        find_conn_cliques(prod, isconn; kwargs...) : find_cliques(prod; kwargs...))
    # delta-Y detection
    mapping = Dict{Edge{T},Edge{T}}()
    for cq in cqstate.cliques
        length(cq) > length(mapping) || continue
        # modprod reverse mapping
        lnnodemap = Dict(div(i - 1, ne(h)) + 1 => mod(i - 1, ne(h)) + 1 for i in cq)  # lg(v) => lh(v)
        emap = Dict(grev[m] => hrev[n] for (m, n) in lnnodemap)  # g(e) => h(e)
        delta_y_correction!(emap, g, h)
        if length(emap) > length(mapping)
            mapping = emap
        end
    end
    return MCSResult(mapping, cqstate.status)
end