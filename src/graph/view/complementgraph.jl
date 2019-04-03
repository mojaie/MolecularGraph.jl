#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

# TODO: deprecated
export
    ComplementGraphView, ComplementDiGraphView,
    complementgraph


struct ComplementGraphView{T<:UndirectedGraph} <: GraphView
    graph::T
end


struct ComplementDiGraphView{T<:DirectedGraph} <: DiGraphView
    graph::T
end

# TODO: use traits
ComplementGraph = Union{ComplementGraphView,ComplementDiGraphView}

complementgraph(graph::UndirectedGraph) = ComplementGraphView(graph)
complementgraph(graph::DirectedGraph) = ComplementDiGraphView(graph)


function getedge(view::ComplementGraphView, idx)
    throw(ErrorException("Edge indexing not supported for this view."))
end

function getedge(view::ComplementGraphView, u, v)
    (v in neighbors(view.graph, u)) || throw(KeyError(v))
    return edgetype(view.graph)(u, v)
end

hasedge(view::ComplementGraphView, u, v) = !hasedge(view.graph, u, v)


function edgesiter(view::ComplementGraphView)
    throw(ErrorException("Edge indexing not supported for this view."))
end

function neighbors(view::ComplementGraphView, idx)
    throw(ErrorException("Edge indexing not supported for this view."))
end

function adjacencies(view::ComplementGraphView, idx)
    ns = nodeset(view.graph)
    pop!(ns, idx)
    return setdiff(ns, adjacencies(view.graph, idx))
end

function incidences(view::ComplementGraphView, idx)
    throw(ErrorException("Edge indexing not supported for this view."))
end

function incidentedges(view::ComplementGraphView, idx)
    throw(ErrorException("Edge indexing not supported for this view."))
end

function edgecount(view::ComplementGraphView)
    n = nodecount(view.graph)
    return n * (n - 1) / 2 - edgecount(view.graph)
end
