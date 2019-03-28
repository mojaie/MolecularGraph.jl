#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    SubgraphView, DiSubgraphView,
    nodesubgraph, edgesubgraph


struct SubgraphView{T<:UndirectedGraph} <: GraphView
    graph::T
    nodes::Set{Int}
    edges::Set{Int}
end


struct DiSubgraphView{T<:DirectedGraph} <: DiGraphView
    graph::T
    nodes::Set{Int}
    edges::Set{Int}
end


# TODO: use traits
Subgraph = Union{SubgraphView,DiSubgraphView}


Base.getindex(view::Subgraph, sym::Symbol) = eval(Expr(:call, sym, view.graph))
Base.getindex(
    graph::Subgraph, k1::Symbol, k2::Symbol, K::Symbol...
) = hcat(eval(Expr(:call, sym, view.graph)) for k in [k1, k2, K...])


"""
    nodesubgraph(graph::UndirectedGraph, nodes) -> SubgraphView
    nodesubgraph(graph::DirectedGraph, nodes) -> DiSubgraphView

Generate node-induced subgraph view.
"""
nodesubgraph(
    graph::UndirectedGraph, nodes
) = SubgraphView(_nodesubgraph(graph, nodes)...)

nodesubgraph(
    graph::DirectedGraph, nodes
) = DiSubgraphView(_nodesubgraph(graph, nodes)...)

function _nodesubgraph(graph, nodes)
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
    edgesubgraph(graph::UndirectedGraph, edges) -> SubgraphView
    edgesubgraph(graph::DirectedGraph, edges) -> DiSubgraphView

Generate edge-induced subgraph view.
"""
function edgesubgraph(graph::UndirectedGraph, edges)
    nodes = Set{Int}()
    for e in edges
        eg = getedge(graph, e)
        push!(nodes, eg.u, eg.v)
    end
    return SubgraphView(graph, nodes, Set(edges))
end

function edgesubgraph(graph::DirectedGraph, edges)
    nodes = Set{Int}()
    for e in edges
        eg = getedge(graph, e)
        push!(nodes, eg.source, eg.target)
    end
    return DiSubgraphView(graph, nodes, Set(edges))
end


function getnode(view::Subgraph, idx)
    (idx in view.nodes) || throw(KeyError(idx))
    return getnode(view.graph, idx)
end


function getedge(view::Subgraph, idx)
    (idx in view.edges) || throw(KeyError(idx))
    return getedge(view.graph, idx)
end

function getedge(view::Subgraph, u, v)
    (u in view.nodes) || throw(KeyError(u))
    (v in view.nodes) || throw(KeyError(v))
    return getedge(view.graph, u, v)
end


nodesiter(view::Subgraph) = [i => getnode(view.graph, i) for i in view.nodes]
edgesiter(view::Subgraph) = [e => getedge(view.graph, e) for e in view.edges]

nodekeys(view::Subgraph) = collect(view.nodes)
edgekeys(view::Subgraph) = collect(view.edges)

nodeset(view::Subgraph) = copy(view.nodes)
edgeset(view::Subgraph) = copy(view.edges)


function neighbors(view::Subgraph, idx)
    (idx in view.nodes) || throw(KeyError(idx))
    return Dict(
        n for n in neighbors(view.graph, idx)
        if n.first in view.nodes && n.second in view.edges
    )
end


function outneighbors(view::DiSubgraphView, idx)
    (idx in view.nodes) || throw(KeyError(idx))
    return Dict(
        n for n in outneighbors(view.graph, idx)
        if n.first in view.nodes && n.second in view.edges
    )
end


function inneighbors(view::DiSubgraphView, idx)
    (idx in view.nodes) || throw(KeyError(idx))
    return Dict(
        n for n in inneighbors(view.graph, idx)
        if n.first in view.nodes && n.second in view.edges
    )
end


nodecount(view::Subgraph) = length(view.nodes)
edgecount(view::Subgraph) = length(view.edges)
