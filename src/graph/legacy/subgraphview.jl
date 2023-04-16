#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    SubgraphView, DiSubgraphView,
    subgraphview, supergraph, nodesubgraph, edgesubgraph


struct SubgraphView{T<:UndirectedGraph} <: UndirectedGraph
    graph::T
    nodes::Set{Int}
    edges::Set{Int}
end
SubgraphView(graph::UndirectedGraph, nodes, edges) = SubgraphView(graph, Set{Int}(nodes), Set{Int}(edges))

Base.:(==)(g1::SubgraphView, g2::SubgraphView) = g1.graph == g2.graph && g1.nodes == g2.nodes && g1.edges == g2.edges

struct DiSubgraphView{T<:DirectedGraph} <: DirectedGraph
    graph::T
    nodes::Set{Int}
    edges::Set{Int}
end
DiSubgraphView(graph::UndirectedGraph, nodes, edges) = DiSubgraphView(graph, Set{Int}(nodes), Set{Int}(edges))

Base.:(==)(g1::DiSubgraphView, g2::DiSubgraphView) = g1.graph == g2.graph && g1.nodes == g2.nodes && g1.edges == g2.edges


"""
    subgraphview(graph::UndirectedGraph, nodes, edges) -> SubgraphView
    subgraphview(graph::DirectedGraph, nodes, edges) -> DiSubgraphView

Generate subgraph view with the given nodes and edges.
"""
subgraphview(graph::UndirectedGraph, nodes, edges) = SubgraphView(graph, nodes, edges)
subgraphview(graph::DirectedGraph, nodes, edges) = DiSubgraphView(graph, nodes, edges)


"""
    supergraph(view::SubgraphView) -> UndirectedGraph
    supergraph(view::DiSubgraphView) -> DirectedGraph

Returns the graph on which the view is based.
"""
supergraph(view::SubgraphView) = view.graph
supergraph(view::DiSubgraphView) = view.graph


function nodesubgraphedges(graph::AbstractGraph, nodes)
    edges = Set{Int}()
    for n in nodes
        for (ninc, nadj) in neighbors(graph, n)
            nadj in nodes && push!(edges, ninc)
        end
    end
    return edges
end


"""
    nodesubgraph(graph::UndirectedGraph, nodes) -> SubgraphView
    nodesubgraph(graph::DirectedGraph, nodes) -> DiSubgraphView

Generate node-induced subgraph view.
"""
function nodesubgraph(graph, nodes)
    edges = nodesubgraphedges(graph, nodes)
    return subgraphview(graph, Set(nodes), edges)
end


function edgesubgraphnodes(graph::AbstractGraph, edges)
    nodes = Set{Int}()
    for e in edges
        union!(nodes, getedge(graph, e))
    end
    return nodes
end


"""
    edgesubgraph(graph::UndirectedGraph, edges) -> SubgraphView
    edgesubgraph(graph::DirectedGraph, edges) -> DiSubgraphView

Generate edge-induced subgraph view.
"""
function edgesubgraph(graph, edges)
    nodes = edgesubgraphnodes(graph, edges)
    return subgraphview(graph, nodes, Set(edges))
end


function neighbors(view::Union{SubgraphView,DiSubgraphView}, idx)
    (idx in view.nodes) || throw(KeyError(idx))
    return Dict(
        n for n in neighbors(view.graph, idx)
        if n.first in view.edges && n.second in view.nodes
    )
end

function outneighbors(view::DiSubgraphView, idx)
    (idx in view.nodes) || throw(KeyError(idx))
    return Dict(
        n for n in outneighbors(view.graph, idx)
        if n.first in view.edges && n.second in view.nodes
    )
end

function inneighbors(view::DiSubgraphView, idx)
    (idx in view.nodes) || throw(KeyError(idx))
    return Dict(
        n for n in inneighbors(view.graph, idx)
        if n.first in view.edges && n.second in view.nodes
    )
end


function getedge(view::Union{SubgraphView,DiSubgraphView}, idx)
    idx in view.edges || throw(KeyError(idx))
    return getedge(view.graph, idx)
end

function findedgekey(view::Union{SubgraphView,DiSubgraphView}, u, v)
    u in view.nodes || return
    v in view.nodes || return
    return findedgekey(view.graph, u, v)
end

function nodeattr(view::Union{SubgraphView,DiSubgraphView}, idx)
    idx in view.nodes || throw(KeyError(idx))
    return nodeattr(view.graph, idx)
end

function edgeattr(view::Union{SubgraphView,DiSubgraphView}, idx)
    idx in view.edges || throw(KeyError(idx))
    return edgeattr(view.graph, idx)
end


nodeset(view::Union{SubgraphView,DiSubgraphView}) = copy(view.nodes)
edgeset(view::Union{SubgraphView,DiSubgraphView}) = copy(view.edges)

nodecount(view::Union{SubgraphView,DiSubgraphView}) = length(view.nodes)
edgecount(view::Union{SubgraphView,DiSubgraphView}) = length(view.edges)

nodeattrs(view::Union{SubgraphView,DiSubgraphView}) = nodeattrs(view.graph)
edgeattrs(view::Union{SubgraphView,DiSubgraphView}) = edgeattrs(view.graph)

nodeattrtype(view::Union{SubgraphView,DiSubgraphView}) = nodeattrtype(view.graph)
edgeattrtype(view::Union{SubgraphView,DiSubgraphView}) = edgeattrtype(view.graph)


hascachestore(view::Union{SubgraphView,DiSubgraphView}) = hascachestore(view.graph)
hascache(view::Union{SubgraphView,DiSubgraphView}, key) = hascache(view.graph, key)
getcache(view::Union{SubgraphView,DiSubgraphView}, key) = getcache(view.graph, key)
setcache!(view::Union{SubgraphView,DiSubgraphView}, key; kwargs...) = setcache!(view.graph, key; kwargs...)
clearcache!(view::Union{SubgraphView,DiSubgraphView}) = clearcache!(view.graph)