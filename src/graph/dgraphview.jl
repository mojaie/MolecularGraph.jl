#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    GDGraphView,
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
    similarmap,
    ReverseGraph,
    edgesiter,
    successors,
    predecessors


struct GDGraphView{G<:AbstractDGraph} <: DGraphView
    graph::G
    nodes::Set{Int}
    edges::Set{Int}
end


function nodesubgraph(graph::AbstractDGraph, nodes)
    """ Node induced subgraph"""
    edges = Set{Int}()
    for n in nodes
        for (nbr, a) in successors(graph, n)
            if nbr in nodes
                push!(edges, a)
            end
        end
    end
    return GDGraphView(graph, Set(nodes), edges)
end


function edgesubgraph(graph::AbstractDGraph, edges)
    """ Edge induced subgraph"""
    nodes = Set{Int}()
    for a in edges
        ag = getedge(graph, a)
        push!(nodes, ag.source, ag.target)
    end
    return GDGraphView(graph, nodes, Set(edges))
end


function getnode(view::DGraphView, idx)
    if !(idx in view.nodes)
        throw(OperationError("Missing node: $(idx)"))
    end
    getnode(view.graph, idx)
end


function getedge(view::DGraphView, idx)
    if !(idx in view.edges)
        throw(OperationError("Missing edge: $(idx)"))
    end
    getedge(view.graph, idx)
end

function getedge(view::DGraphView, source, target)
    if !(source in view.nodes)
        throw(OperationError("Missing node: $(source)"))
    elseif !(target in view.nodes)
        throw(OperationError("Missing node: $(target)"))
    end
    getedge(view.graph, source, target)
end


nodesiter(view::DGraphView) = [i => getnode(view.graph, i) for i in view.nodes]

edgesiter(view::DGraphView) = [
    a => getedge(view.graph, a) for a in view.edges]

nodekeys(view::DGraphView) = copy(view.nodes)

edgekeys(view::DGraphView) = copy(view.edges)


function successors(view::DGraphView, idx)
    if !(idx in view.nodes)
        throw(OperationError("Missing node: $(idx)"))
    end
    return Dict(
        n for n in successors(view.graph, idx)
        if n.first in view.nodes && n.second in view.edges
    )
end

function predecessors(view::DGraphView, idx)
    if !(idx in view.nodes)
        throw(OperationError("Missing node: $(idx)"))
    end
    return Dict(
        n for n in predecessors(view.graph, idx)
        if n.first in view.nodes && n.second in view.edges
    )
end

nodecount(view::DGraphView) = length(view.nodes)

edgecount(view::DGraphView) = length(view.edges)

similarmap(view::DGraphView) = similarmap(view.graph)



struct ReverseGraph{G<:AbstractDGraph} <: DGraphView
    graph::G
end

function getedge(view::ReverseGraph, idx)
    e = getedge(view.graph, idx)
    return connect(e, e.target, e.source)
end

getedge(view::ReverseGraph, s, t) = connect(getedge(view.graph, s, t), t, s)

edgesiter(view::ReverseGraph) = [
    i => connect(e, e.target, e.source) for (i, e) in edgesiter(view.graph)]

successors(view::ReverseGraph, idx) = view.graph.predecessors[idx]

predecessors(view::ReverseGraph, idx) = view.graph.successors[idx]
