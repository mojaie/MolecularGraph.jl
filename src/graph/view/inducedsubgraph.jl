#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    UDSubgraph, DSubgraph, SubgraphView,
    nodesubgraph, edgesubgraph,
    getnode, getedge,
    nodesiter, edgesiter,
    nodekeys, edgekeys,
    neighbors, successors, predecessors,
    nodecount, edgecount


struct UDSubgraph{T<:UDGraph} <: UndirectedGraphView
    graph::T
    nodes::Set{Int}
    edges::Set{Int}
end


struct DSubgraph{T<:DGraph} <: DirectedGraphView
    graph::T
    nodes::Set{Int}
    edges::Set{Int}
end


# TODO: use traits
SubgraphView = Union{UDSubgraph,DSubgraph}


"""
    nodesubgraph(graph::AbstractGraph, nodes)

Generate node-induced subgraph view.
"""
nodesubgraph(g::UDGraph, nodes) = UDSubgraph(_nodesubgraph(g, nodes)...)
nodesubgraph(g::DGraph, nodes) = DSubgraph(_nodesubgraph(g, nodes)...)
function _nodesubgraph(graph, nodes)
    """ Node induced subgraph"""
    edges = Set{Int}()
    for n in nodes
        for (nbr, e) in neighbors(graph, n)
            if nbr in nodes
                push!(edges, e)
            end
        end
    end
    return (graph, Set(nodes), edges)
end


"""
    edgesubgraph(graph::AbstractGraph, edges)

Generate edge-induced subgraph view.
"""
function edgesubgraph(graph::UDGraph, edges)
    """ Edge induced subgraph"""
    nodes = Set{Int}()
    for e in edges
        eg = getedge(graph, e)
        push!(nodes, eg.u, eg.v)
    end
    return UDSubgraph(graph, nodes, Set(edges))
end

function edgesubgraph(graph::DGraph, edges)
    """ Edge induced subgraph"""
    nodes = Set{Int}()
    for e in edges
        eg = getedge(graph, e)
        push!(nodes, eg.source, eg.target)
    end
    return DSubgraph(graph, nodes, Set(edges))
end


function getnode(view::SubgraphView, idx)
    (idx in view.nodes) || throw(KeyError(idx))
    getnode(view.graph, idx)
end


function getedge(view::SubgraphView, idx)
    (idx in view.edges) || throw(KeyError(idx))
    getedge(view.graph, idx)
end

function getedge(view::SubgraphView, u, v)
    (u in view.nodes) || throw(KeyError(u))
    (v in view.nodes) || throw(KeyError(v))
    getedge(view.graph, u, v)
end


nodesiter(view::SubgraphView) = (
    i => getnode(view.graph, i) for i in view.nodes)
edgesiter(view::SubgraphView) = (
    e => getedge(view.graph, e) for e in view.edges)

nodekeys(view::SubgraphView) = copy(view.nodes)
edgekeys(view::SubgraphView) = copy(view.edges)


function neighbors(view::SubgraphView, idx)
    (idx in view.nodes) || throw(KeyError(idx))
    return Dict(
        n for n in neighbors(view.graph, idx)
        if n.first in view.nodes && n.second in view.edges
    )
end


function successors(view::DSubgraph, idx)
    (idx in view.nodes) || throw(KeyError(idx))
    return Dict(
        n for n in successors(view.graph, idx)
        if n.first in view.nodes && n.second in view.edges
    )
end


function predecessors(view::DSubgraph, idx)
    (idx in view.nodes) || throw(KeyError(idx))
    return Dict(
        n for n in predecessors(view.graph, idx)
        if n.first in view.nodes && n.second in view.edges
    )
end


nodecount(view::SubgraphView) = length(view.nodes)
edgecount(view::SubgraphView) = length(view.edges)
