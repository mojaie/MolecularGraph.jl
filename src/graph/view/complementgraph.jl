#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    UDComplementGraph, DComplementGraph,
    ComplementGraphView,
    edgecount


struct UDComplementGraph{T<:UDGraph} <: UndirectedGraphView
    graph::T
end


struct DComplementGraph{T<:DGraph} <: DirectedGraphView
    graph::T
end

# TODO: use traits
ComplementGraphView = Union{UDComplementGraph,DComplementGraph}


function getedge(view::UDComplementGraph, idx)
    throw(ErrorException("Edge indexing not supported for this view."))
end

function getedge(view::UDComplementGraph, u, v)
    # TODO: hasedge
    (v in neighbors(view.graph, u)) || throw(KeyError(v))
    e = getedge(view.graph, u, v)
    return similaredge(e, u, v)
end

function edgesiter(view::UDComplementGraph)
    throw(ErrorException("Edge indexing not supported for this view."))
end

function neighbors(view::UDComplementGraph, idx)
    throw(ErrorException("Edge indexing not supported for this view."))
end

function neighborset(view::UDComplementGraph, idx)
    ns = nodeset(view.graph)
    pop!(ns, idx)
    return setdiff(ns, neighborset(view.graph, idx))
end

function neighboredgeset(view::UDComplementGraph, idx)
    throw(ErrorException("Edge indexing not supported for this view."))
end

function neighboredges(view::UDComplementGraph, idx)
    throw(ErrorException("Edge indexing not supported for this view."))
end

function edgecount(view::UDComplementGraph)
    n = nodecount(view.graph)
    return n * (n - 1) / 2 - edgecount(view.graph)
end
