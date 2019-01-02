#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    UDSubgraph,
    DSubgraph,
    SubgraphView,
    nodesubgraph,
    edgesubgraph,
    getnode,
    getedge,
    nodesiter,
    edgesiter,
    nodekeys,
    edgekeys,
    neighbors,
    nodecount,
    edgecount,
    similarmap


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


edgesubgraph(g::UDGraph, nodes) = UDSubgraph(_edgesubgraph(g, nodes)...)
edgesubgraph(g::DGraph, nodes) = DSubgraph(_edgesubgraph(g, nodes)...)
function _edgesubgraph(graph, edges)
    """ Edge induced subgraph"""
    nodes = Set{Int}()
    for e in edges
        eg = getedge(graph, e)
        push!(nodes, eg.u, eg.v)
    end
    return (graph, nodes, Set(edges))
end


function getnode(view::SubgraphView, idx)
    (idx in view.nodes) || throw(OperationError("Missing node: $(idx)"))
    getnode(view.graph, idx)
end


function getedge(view::SubgraphView, idx)
    (idx in view.edges) || throw(OperationError("Missing edge: $(idx)"))
    getedge(view.graph, idx)
end

function getedge(view::SubgraphView, u, v)
    (u in view.nodes) || throw(OperationError("Missing node: $(u)"))
    (v in view.nodes) || throw(OperationError("Missing node: $(v)"))
    getedge(view.graph, u, v)
end


nodesiter(view::SubgraphView) = (
    i => getnode(view.graph, i) for i in view.nodes)
edgesiter(view::SubgraphView) = (
    e => getedge(view.graph, e) for e in view.edges)

nodekeys(view::SubgraphView) = copy(view.nodes)
edgekeys(view::SubgraphView) = copy(view.edges)


function neighbors(view::SubgraphView, idx)
    (idx in view.nodes) || throw(OperationError("Missing node: $(idx)"))
    return Dict(
        n for n in neighbors(view.graph, idx)
        if n.first in view.nodes && n.second in view.edges
    )
end


function successors(view::DSubgraph, idx)
    (idx in view.nodes) || throw(OperationError("Missing node: $(idx)"))
    return Dict(
        n for n in successors(view.graph, idx)
        if n.first in view.nodes && n.second in view.edges
    )
end

function predecessors(view::DSubgraph, idx)
    (idx in view.nodes) || throw(OperationError("Missing node: $(idx)"))
    return Dict(
        n for n in predecessors(view.graph, idx)
        if n.first in view.nodes && n.second in view.edges
    )
end


nodecount(view::SubgraphView) = length(view.nodes)
edgecount(view::SubgraphView) = length(view.edges)

similarmap(view::SubgraphView) = similarmap(view.graph)


struct ReverseGraph{T<:DGraph} <: DirectedGraphView
    graph::T
end

function getedge(view::ReverseGraph, idx)
    e = getedge(view.graph, idx)
    return connect(e, e.target, e.source)
end

getedge(view::ReverseGraph, s, t) = connect(getedge(view.graph, s, t), t, s)

edgesiter(view::ReverseGraph) = (
    i => connect(e, e.target, e.source) for (i, e) in edgesiter(view.graph))

successors(view::ReverseGraph, idx) = view.graph.predecessors[idx]

predecessors(view::ReverseGraph, idx) = view.graph.successors[idx]
