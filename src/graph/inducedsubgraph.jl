#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    SubgraphView, DiSubgraphView,
    supergraph, nodesubgraph, edgesubgraph


struct SubgraphView{T<:UndirectedGraph} <: UndirectedGraph
    graph::T
    nodes::Set{Int}
    edges::Set{Int}
end

struct DiSubgraphView{T<:DirectedGraph} <: DirectedGraph
    graph::T
    nodes::Set{Int}
    edges::Set{Int}
end


# TODO: use trait
AllSubgraphView = Union{SubgraphView,DiSubgraphView}


supergraph(view::AllSubgraphView) = view.graph


"""
    nodesubgraph(graph::UndirectedGraph, nodes::Set{Int}) -> SubgraphView
    nodesubgraph(graph::DirectedGraph, nodes::Set{Int}) -> DiSubgraphView

Generate node-induced subgraph view.
"""
function nodesubgraph(graph::AbstractGraph, nodes::Set{Int})
    # TODO: test for updating in nodes
    edges = Set{Int}()
    for n in nodes
        for (ninc, nadj) in neighbors(graph, n)
            if nadj in nodes
                push!(edges, ninc)
            end
        end
    end
    if graph isa UndirectedGraph
        return SubgraphView(graph, Set(nodes), edges)
    elseif graph isa DirectedGraph
        return DiSubgraphView(graph, Set(nodes), edges)
    end
end


"""
    edgesubgraph(graph::UndirectedGraph, edges::Set{Int}) -> SubgraphView
    edgesubgraph(graph::DirectedGraph, edges::Set{Int}) -> DiSubgraphView

Generate edge-induced subgraph view.
"""
function edgesubgraph(graph::AbstractGraph, edges::Set{Int})
    # TODO: test for updating in edges
    nodes = Set{Int}()
    for e in edges
        union!(nodes, getedge(graph, e))
    end
    if graph isa UndirectedGraph
        return SubgraphView(graph, nodes, Set(edges))
    elseif graph isa DirectedGraph
        return DiSubgraphView(graph, nodes, Set(edges))
    end
end



function neighbors(view::AllSubgraphView, idx::Int)
    (idx in view.nodes) || throw(KeyError(idx))
    return Dict(
        n for n in neighbors(view.graph, idx)
        if n.first in view.edges && n.second in view.nodes
    )
end

function outneighbors(view::DiSubgraphView, idx::Int)
    (idx in view.nodes) || throw(KeyError(idx))
    return Dict(
        n for n in outneighbors(view.graph, idx)
        if n.first in view.edges && n.second in view.nodes
    )
end

function inneighbors(view::DiSubgraphView, idx::Int)
    (idx in view.nodes) || throw(KeyError(idx))
    return Dict(
        n for n in inneighbors(view.graph, idx)
        if n.first in view.edges && n.second in view.nodes
    )
end


function getedge(view::AllSubgraphView, idx::Int)
    idx in view.edges || throw(KeyError(idx))
    return getedge(view.graph, idx)
end

function hasedge(view::AllSubgraphView, u::Int, v::Int)
    u in view.nodes || return false
    v in view.nodes || return false
    return findedgekey(view.graph, u, v) !== nothing
end

function nodeattr(view::AllSubgraphView, idx::Int)
    idx in view.nodes || throw(KeyError(idx))
    return nodeattr(view.graph, idx)
end

function edgeattr(view::AllSubgraphView, idx::Int)
    idx in view.edges || throw(KeyError(idx))
    return edgeattr(view.graph, idx)
end

function edgeattr(view::AllSubgraphView, u::Int, v::Int)
    u in view.nodes || return nothing
    v in view.nodes || return nothing
    k = findedgekey(view.graph, u, v)
    return k === nothing ? nothing : edgeattr(view.graph, k)
end


nodeset(view::AllSubgraphView) = copy(view.nodes)
edgeset(view::AllSubgraphView) = copy(view.edges)

nodecount(view::Union{SubgraphView,DiSubgraphView}) = length(view.nodes)
edgecount(view::Union{SubgraphView,DiSubgraphView}) = length(view.edges)

nodeattrs(view::AllSubgraphView) = nodeattrs(view.graph)
edgeattrs(view::AllSubgraphView) = edgeattrs(view.graph)

nodeattrtype(view::AllSubgraphView) = nodeattrtype(view.graph)
edgeattrtype(view::AllSubgraphView) = edgeattrtype(view.graph)
