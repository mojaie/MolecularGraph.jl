#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    ReverseGraph,
    getedge,
    edgesiter,
    successors, predecessors


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
