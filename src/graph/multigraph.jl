#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    MultiUDGraph,
    getnode, getedge, hasedge,
    nodesiter, edgesiter,
    nodekeys, edgekeys,
    neighbors,
    updatenode!, updateedge!,
    unlinknode!, unlinkedge!


struct MultiUDGraph{N<:AbstractNode,E<:UndirectedEdge} <: UndirectedGraph
    nodes::Dict{Int,N}
    edges::Dict{Int,E}
    adjacency::Dict{Int,Dict{Int,Set{Int}}}

    function MultiUDGraph{N,E}() where {N<:AbstractNode,E<:UndirectedEdge}
        new(Dict(), Dict(), Dict())
    end
end

function MultiUDGraph(nodes, edges)
    graph = MultiUDGraph{Node,Edge}()
    for node in nodes
        updatenode!(graph, Node(), node)
    end
    for (i, edge) in enumerate(edges)
        updateedge!(graph, Edge(edge...), i)
    end
    return graph
end

getnode(graph::MultiUDGraph, idx) = graph.nodes[idx]

getedge(graph::MultiUDGraph, idx) = graph.edges[idx]

hasedge(graph::MultiUDGraph, u, v) = haskey(graph.adjacency[u], v)

nodesiter(graph::MultiUDGraph) = graph.nodes

nodekeys(graph::MultiUDGraph) = Set(keys(graph.nodes))

edgesiter(graph::MultiUDGraph) = graph.edges

edgekeys(graph::MultiUDGraph) = Set(keys(graph.edges))

function neighbors(graph::MultiUDGraph, idx)
    nbrs = Pair{Int,Int}[]
    for (n, edges) in graph.adjacency[idx]
        for e in edges
            push!(nbrs, n => e)
        end
    end
    return nbrs
end

neighboredgekeys(g::MultiUDGraph, i) = Set([nbr.second for nbr in neighbors(g, i)])


function updatenode!(graph::MultiUDGraph, node, idx)
    """Add or update a node"""
    graph.nodes[idx] = node
    if !(idx in keys(graph.adjacency))
        graph.adjacency[idx] = Dict()
    end
    return
end


function updateedge!(graph::MultiUDGraph, edge, idx)
    """Add or update an edge"""
    nodes = nodekeys(graph)
    (edge.u in nodes) || throw(KeyError(edge.u))
    (edge.v in nodes) || throw(KeyError(edge.v))
    graph.edges[idx] = edge
    if !haskey(graph.adjacency[edge.u], edge.v)
        graph.adjacency[edge.u][edge.v] = Set{Int}()
    end
    push!(graph.adjacency[edge.u][edge.v], idx)
    if !haskey(graph.adjacency[edge.v], edge.u)
        graph.adjacency[edge.v][edge.u] = Set{Int}()
    end
    push!(graph.adjacency[edge.v][edge.u], idx)
    return
end

function updateedge!(G::MultiUDGraph, edge, u, v)
    throw(ErrorException(
        "The index of the edge to be updated should be specified."))
end


function unlinknode!(graph::MultiUDGraph, idx)
    """Remove a node and its connecting edges"""
    (idx in nodekeys(graph)) || throw(KeyError(idx))
    for (n, edges) in neighbors(graph, idx)
        for e in edges
            delete!(graph.edges, e)
        end
        delete!(graph.adjacency[n], idx)
    end
    delete!(graph.nodes, idx)
    delete!(graph.adjacency, idx)
    return
end


function unlinkedge!(graph::MultiUDGraph, u, v)
    throw(ErrorException(
        "The index of the edge to be deleted should be specified."))
end

function unlinkedge!(graph::MultiUDGraph, idx)
    """Remove an edge"""
    (idx in keys(graph.edges)) || throw(KeyError(idx))
    e = getedge(graph, idx)
    delete!(graph.edges, idx)
    delete!(graph.adjacency[e.u][e.v], idx)
    delete!(graph.adjacency[e.v][e.u], idx)
    isempty(graph.adjacency[e.u]) && delete!(graph.adjacency, e.u)
    isempty(graph.adjacency[e.v]) && delete!(graph.adjacency, e.v)
    return
end
