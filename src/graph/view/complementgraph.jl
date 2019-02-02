#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    UDComplementGraph, DComplementGraph,
    ComplementGraphView,
    getedge,
    edgesiter,
    edgekeys,
    neighbors,
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
    return connect(e, u, v)
end

function edgesiter(view::UDComplementGraph)
    throw(ErrorException("Edge indexing not supported for this view."))
end

function neighbors(view::UDComplementGraph, idx)
    throw(ErrorException("Edge indexing not supported for this view."))
end

function neighborkeys(view::UDComplementGraph, idx)
    ns = nodekeys(view.graph)
    pop!(ns, idx)
    return setdiff(ns, neighborkeys(view.graph, idx))
end

function neighboredgekeys(view::UDComplementGraph, idx)
    throw(ErrorException("Edge indexing not supported for this view."))
end

function neighboredges(view::UDComplementGraph, idx)
    throw(ErrorException("Edge indexing not supported for this view."))
end

function edgecount(view::UDComplementGraph)
    n = nodecount(view.graph)
    return n * (n - 1) / 2 - edgecount(view.graph)
end
