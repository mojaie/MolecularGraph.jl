#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    ReverseGraphView, reversegraph


struct ReverseGraphView{T<:DirectedGraph} <: DiGraphView
    graph::T
end

reversegraph(graph::DirectedGraph) = ReverseGraphView(graph)


function getedge(view::ReverseGraphView, idx)
    edge = getedge(view.graph, idx)
    return reverseedge(edge)
end

getedge(
    view::ReverseGraphView, s, t
) = getedge(view, outneighbors(view.graph, s)[t])

edgesiter(
    view::ReverseGraphView
) = [i => reverseedge(e) for (i, e) in edgesiter(view.graph)]

function reverseedge(edge)
    return setnodes(edge, edge.target, edge.source)
end

outneighbors(view::ReverseGraphView, idx) = view.graph.inneighbors[idx]
inneighbors(view::ReverseGraphView, idx) = view.graph.outneighbors[idx]
