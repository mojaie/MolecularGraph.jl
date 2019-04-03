#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

# TODO: deprecated

export
    PlainMultiGraph, multigraph


struct PlainMultiGraph <: AbstractGraph
    neighbormap::Vector{Dict{Int,Set{Int}}}
    edges::Vector{Edge}
end


"""
    plainmultigraph() -> PlainMultiGraph()

Generate empty `MultiGraph`.
"""
multigraph() = MultiGraph([], [], [])


"""
    multigraph(size::Int, edges) -> MultiGraph

Generate `MultiGraph` that have given nodes and edges represented by the
list of node indices in integer and the list of pairs of node indices,
respectively.
"""
function multigraph(size::Int, edges)
    graph = MultiGraph()
    append!(graph.nodes, [PlainNode() for i in 1:size])
    append!(graph.neighbormap, [Dict() for i in 1:size])
    for (i, (u, v)) in enumerate(edges)
        push!(graph.edges, Edge(u, v))
        haskey(graph.neighbormap[u], v) || (graph.neighbormap[u][v] = Set())
        haskey(graph.neighbormap[v], u) || (graph.neighbormap[v][u] = Set())
        push!(graph.neighbormap[u][v], i)
        push!(graph.neighbormap[v][u], i)
    end
    return graph
end


getnode(graph::MultiGraph, idx) = graph.nodes[idx]
getedge(graph::MultiGraph, idx) = graph.edges[idx]
hasedge(graph::MultiGraph, u, v) = haskey(graph.neighbormap[u], v)


nodekeys(graph::MultiGraph) = collect(1:nodecount(graph))
edgekeys(graph::MultiGraph) = collect(1:edgecount(graph))
nodevalues(graph::MultiGraph) = graph.nodes
edgevalues(graph::MultiGraph) = graph.edges
nodesiter(graph::MultiGraph) = enumerate(graph.nodes)
edgesiter(graph::MultiGraph) = enumerate(graph.edges)
nodeset(graph::MultiGraph) = Set(1:nodecount(graph))
edgeset(graph::MultiGraph) = Set(1:edgecount(graph))

nodecount(graph::MultiGraph) = length(graph.nodes)
edgecount(graph::MultiGraph) = length(graph.edges)


function neighbors(graph::MultiGraph, idx)
    nbrs = Pair{Int,Int}[]
    for (n, edges) in graph.neighbormap[idx]
        for e in edges
            push!(nbrs, n => e)
        end
    end
    return nbrs
end

adjacencies(g::MultiGraph, i) = Set([nbr.first for nbr in neighbors(g, i)])
incidences(g::MultiGraph, i) = Set([nbr.second for nbr in neighbors(g, i)])


function updatenode!(graph::MultiGraph, node, idx)
    idx > nodecount(graph) + 1 && throw(DomainError(idx))
    if idx == nodecount(graph) + 1
        updatenode!(graph, node)
        return
    end
    graph.nodes[idx] = node
    return
end

function updatenode!(graph::MultiGraph, node)
    push!(graph.nodes, node)
    push!(graph.neighbormap, Dict())
    return
end


function updateedge!(graph::MultiGraph, edge, idx)
    idx > edgecount(graph) + 1 && throw(DomainError(idx))
    if idx == edgecount(graph) + 1
        updateedge!(graph, edge)
        return
    end
    ncnt = nodecount(graph)
    edge.u > ncnt && throw(DomainError(edge.u))
    edge.v > ncnt && throw(DomainError(edge.v))
    old = getedge(graph, idx)
    delete!(graph.neighbormap[old.u][old.v], idx)
    delete!(graph.neighbormap[old.v][old.u], idx)
    if isempty(graph.neighbormap[old.u][old.v])
        delete!(graph.neighbormap[old.u], old.v)
    end
    if isempty(graph.neighbormap[old.v][old.u])
        delete!(graph.neighbormap[old.v], old.u)
    end
    graph.edges[idx] = edge
    haskey(graph.neighbormap[u], v) || (graph.neighbormap[u][v] = Set{Int}())
    haskey(graph.neighbormap[v], u) || (graph.neighbormap[v][u] = Set{Int}())
    push!(graph.neighbormap[edge.u][edge.v], idx)
    push!(graph.neighbormap[edge.v][edge.u], idx)
    return
end

function updateedge!(graph::MultiGraph, edge)
    ncnt = nodecount(graph)
    edge.u > ncnt && throw(DomainError(edge.u))
    edge.v > ncnt && throw(DomainError(edge.v))
    i = edgecount(graph) + 1
    push!(graph.edges, edge)
    haskey(graph.neighbormap[u], v) || (graph.neighbormap[u][v] = Set{Int}())
    haskey(graph.neighbormap[v], u) || (graph.neighbormap[v][u] = Set{Int}())
    push!(graph.neighbormap[edge.u][edge.v], i)
    push!(graph.neighbormap[edge.v][edge.u], i)
    return
end

function updateedge!(G::MultiGraph, edge, u, v)
    throw(ErrorException(
        "The index of the edge to be updated should be specified."))
end

# TODO: unlinknodes, unlinkedges

nodetype(graph::MultiGraph) = eltype(graph.nodes)
edgetype(graph::MultiGraph) = eltype(graph.edges)
