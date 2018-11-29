#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    GUGraphView,
    inducedsubgraph,
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


struct GUGraphView{G<:AbstractUGraph} <: UGraphView
    graph::G
    nodes::Set{Int}
    edges::Set{Int}
end


function inducedsubgraph(graph::AbstractUGraph, nodes)
    edges = Set{Int}()
    for n in nodes
        for (nbr, e) in neighbors(graph, n)
            if nbr in nodes
                push!(edges, e)
            end
        end
    end
    return GUGraphView(graph, Set(nodes), edges)
end


function getnode(view::UGraphView, idx)
    if !(idx in view.nodes)
        throw(OperationError("Missing node: $(idx)"))
    end
    getnode(view.graph, idx)
end


function getedge(view::UGraphView, idx)
    if !(idx in view.edges)
        throw(OperationError("Missing edge: $(idx)"))
    end
    getedge(view.graph, idx)
end

function getedge(view::UGraphView, u, v)
    if !(u in view.nodes)
        throw(OperationError("Missing node: $(u)"))
    elseif !(v in view.nodes)
        throw(OperationError("Missing node: $(v)"))
    end
    getedge(view.graph, u, v)
end


nodesiter(view::UGraphView) = [i => getnode(view.graph, i) for i in view.nodes]

edgesiter(view::UGraphView) = [e => getedge(view.graph, e) for e in view.edges]

nodekeys(view::UGraphView) = copy(view.nodes)

edgekeys(view::UGraphView) = copy(view.edges)


function neighbors(view::UGraphView, idx)
    if !(idx in view.nodes)
        throw(OperationError("Missing node: $(idx)"))
    end
    return Dict(
        n for n in neighbors(view.graph, idx)
        if n.first in view.nodes && n.second in view.edges
    )
end


nodecount(view::UGraphView) = length(view.nodes)

edgecount(view::UGraphView) = length(view.edges)

similarmap(view::UGraphView) = similarmap(view.graph)
