#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    MapMultiGraph, mapmultigraph


struct MapMultiGraph{N<:AbstractNode,E<:UndirectedEdge} <: MultiGraph
    nodes::Dict{Int,N}
    edges::Dict{Int,E}
    neighbormap::Dict{Int,Dict{Int,Set{Int}}}

    function MapMultiGraph{N,E}() where {N<:AbstractNode,E<:UndirectedEdge}
        new(Dict(), Dict(), Dict())
    end
end


"""
    mapmultigraph(::Type{N}, ::Type{E}
        ) where {N<:AbstractNode,E<:UndirectedEdge} -> MapMultiGraph{N,E}()

Generate empty `MapMultiGraph` that have nodes and edges with the given types.
"""
mapmultigraph(::Type{N}, ::Type{E}
    ) where {N<:AbstractNode,E<:UndirectedEdge} = MapMultiGraph{N,E}()

"""
    mapmultigraph(nodes, edges) -> MapMultiGraph{Node,Edge}

Generate `MapMultiGraph` that have given nodes and edges represented by the
list of node indices in integer and the list of pairs of node indices,
respectively.
"""
function mapmultigraph(nodes, edges)
    graph = MapMultiGraph{Node,Edge}()
    for node in nodes
        graph.nodes[node] = Node()
        graph.neighbormap[node] = Dict()
    end
    for (i, edge) in enumerate(edges)
        (u, v) = edge
        graph.edges[i] = Edge(u, v)
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

nodesiter(graph::MultiGraph) = graph.nodes
edgesiter(graph::MultiGraph) = graph.edges
nodekeys(graph::MultiGraph) = collect(keys(graph.nodes))
edgekeys(graph::MultiGraph) = collect(keys(graph.edges))
nodeset(graph::MultiGraph) = Set(keys(graph.nodes))
edgeset(graph::MultiGraph) = Set(keys(graph.edges))


function neighbors(graph::MultiGraph, idx)
    nbrs = Pair{Int,Int}[]
    for (n, edges) in graph.neighbormap[idx]
        for e in edges
            push!(nbrs, n => e)
        end
    end
    return nbrs
end

incidences(g::MultiGraph, i) = Set([nbr.second for nbr in neighbors(g, i)])


function updatenode!(graph::MultiGraph, node, idx)
    graph.nodes[idx] = node
    if !haskey(graph.neighbormap, idx)
        graph.neighbormap[idx] = Dict()
    end
    return
end

function updatenode!(graph::MultiGraph, node)
    i = maximum(nodeset(graph)) + 1
    graph.nodes[i] = node
    graph.neighbormap[i] = Dict()
    return
end


function updateedge!(graph::MultiGraph, edge, idx)
    nodes = nodeset(graph)
    (u, v) = (edge.u, edge.v)
    (u in nodes) || throw(KeyError(u))
    (v in nodes) || throw(KeyError(v))
    if haskey(graph.edges, idx)
        old = getedge(graph, idx)
        delete!(graph.neighbormap[old.u][old.v], idx)
        delete!(graph.neighbormap[old.v][old.u], idx)
        if isempty(graph.neighbormap[old.u][old.v])
            delete!(graph.neighbormap[old.u], old.v)
        end
        if isempty(graph.neighbormap[old.v][old.u])
            delete!(graph.neighbormap[old.v], old.u)
        end
    end
    graph.edges[idx] = edge
    if !haskey(graph.neighbormap[u], v)
        graph.neighbormap[u][v] = Set{Int}()
    end
    if !haskey(graph.neighbormap[v], u)
        graph.neighbormap[v][u] = Set{Int}()
    end
    push!(graph.neighbormap[u][v], idx)
    push!(graph.neighbormap[v][u], idx)
    return
end

updateedge!(graph::MultiGraph, edge
    ) = updateedge!(graph, edge, maximum(edgeset(graph)) + 1)

function updateedge!(G::MultiGraph, edge, u, v)
    throw(ErrorException(
        "The index of the edge to be updated should be specified."))
end


function unlinknode!(graph::MultiGraph, idx)
    (idx in nodekeys(graph)) || throw(KeyError(idx))
    for (n, edges) in neighbors(graph, idx)
        for e in edges
            delete!(graph.edges, e)
        end
        delete!(graph.neighbormap[n], idx)
    end
    delete!(graph.nodes, idx)
    delete!(graph.neighbormap, idx)
    return
end


function unlinkedge!(graph::MultiGraph, u, v)
    throw(ErrorException(
        "The index of the edge to be deleted should be specified."))
end

function unlinkedge!(graph::MultiGraph, idx)
    (idx in keys(graph.edges)) || throw(KeyError(idx))
    e = getedge(graph, idx)
    delete!(graph.edges, idx)
    delete!(graph.neighbormap[e.u][e.v], idx)
    delete!(graph.neighbormap[e.v][e.u], idx)
    isempty(graph.neighbormap[e.u]) && delete!(graph.neighbormap, e.u)
    isempty(graph.neighbormap[e.v]) && delete!(graph.neighbormap, e.v)
    return
end
