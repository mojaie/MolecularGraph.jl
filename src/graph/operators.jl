#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

"""
    induced_subgraph_edges(g::SimpleGraph, node_list::AbstractVector)

Return the node-induced subgraph edges.
"""
function induced_subgraph_edges(g::SimpleGraph, node_list::AbstractVector)
    subg, vmap = induced_subgraph(g, node_list)
    return [u_edge(g, vmap[src(e)], vmap[dst(e)]) for e in edges(subg)]
end


"""
    modularproduct(g::SimpleGraph{T}, h::SimpleGraph{T};
        vmatch=(g1,h1)->true,
        edgefilter=(g1,g2,h1,h2)->has_edge(g,g1,g2)==has_edge(g,h1,h2)) where T

Return the modular product `m` of graphs `g` and `h`, and a mapping whether
the edge is connected or not.
mapping g,h nodes to m nods is f(i, j) = (i - 1) * nv(h) + j and the reverse
mapping is f(i) = (div(i - 1, nv(h)) + 1, mod(i - 1, nv(h))) + 1, where i in g and j in h.
"""
function modular_product(g::SimpleGraph{T}, h::SimpleGraph{T};
            vmatch=(g1,h1)->true,
            edgefilter=(g1,g2,h1,h2)->has_edge(g,g1,g2)==has_edge(g,h1,h2)) where T
    m = SimpleGraph(nv(g) * nv(h))
    connected = Dict{Edge{T},Bool}()  # for connected clique detection
    nv(m) == 0 && return m, connected
    id(i, j) = (i - 1) * nv(h) + j
    for (g1, g2) in combinations(nv(g))
        for (h1, h2) in combinations(nv(h))
            edgefilter(g1, g2, h1, h2) || continue
            if vmatch(g1, h1) && vmatch(g2, h2)
                e = Edge{T}(id(g1, h1), id(g2, h2))  # g1 < g2 -> id(g1, h1) < id(g2, h2)
                add_edge!(m, e)  
                connected[e] = has_edge(g, g1, g2)
            end
            if vmatch(g1, h2) && vmatch(g2, h1)
                e = Edge{T}(id(g1, h2), id(g2, h1))  # g1 < g2 -> id(g1, h2) < id(g2, h1)
                add_edge!(m, e)
                connected[e] = has_edge(g, g1, g2)
            end
        end
    end
    return m, connected
end



"""
    line_graph(G::SimpleGraph) -> SimpleGraph, Dict, Dict

Generate line graph, reverse mapping lg(v) -> g(e) and shared nodes mapping lg(e) -> g(v)
"""
function line_graph(g::SimpleGraph{T}) where T
    l = SimpleGraph{T}(ne(g))
    edge_rank = Dict(e => i for (i, e) in enumerate(edges(g)))
    revmap = Dict(i => e for (i, e) in enumerate(edges(g)))
    sharednode = Dict{Edge{T},T}()
    for i in vertices(g)
        degree(g, i) < 2 && continue
        nbrs = neighbors(g, i)
        for (m, n) in combinations(length(nbrs))
            e = u_edge(T,
                edge_rank[u_edge(T, i, nbrs[m])],
                edge_rank[u_edge(T, i, nbrs[n])]
            )
            add_edge!(l, e)
            sharednode[e] = i
        end
    end
    return l, revmap, sharednode
end


"""
    disjoint_union!(g::T, h::T) where {T<:Graphs.AbstractSimpleGraph} -> Nothing

Add all vertices and edges in graph h to graph g.

Vertices in h will be re-indexed using `i -> i + nv(g)` not to overlap with g vertices.
"""
function disjoint_union!(g::T, h::T) where {T<:Graphs.AbstractSimpleGraph}
    gnv = nv(g)
    for i in 1:nv(h)
        push!(g.fadjlist, eltype(T)[])
        if is_directed(g)
            push!(g.badjlist, eltype(T)[])
        end
    end
    for e in edges(h)
        add_edge!(g, src(e) + gnv, dst(e) + gnv)
    end
    return
end


"""
    disjoint_union(g::T, h::T) where {T<:Graphs.AbstractSimpleGraph} -> Nothing

Return disjoint union of the graph. see `disjoint_union!`
"""
function disjoint_union(g::T, h::T) where {T<:Graphs.AbstractSimpleGraph}
    g_ = copy(g)
    disjoint_union!(g_, h)
    return g_
end