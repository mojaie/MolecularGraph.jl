#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    Edge, similaredge, MapUDGraph, VectorUDGraph,
    getnode, getedge, hasedge,
    nodesiter, edgesiter, nodevector, edgevector,
    neighbors,
    updatenode!, updateedge!,
    unlinknode!, unlinkedge!,
    nodetype, edgetype, similargraph, newgraph


mutable struct Edge <: UndirectedEdge
    u::Int
    v::Int
end

similaredge(e::Edge, u, v) = Edge(u, v)


struct MapUDGraph{N<:AbstractNode,E<:UndirectedEdge} <: UndirectedGraph
    nodes::Dict{Int,N}
    edges::Dict{Int,E}
    adjacency::Dict{Int,Dict{Int,Int}}

    function MapUDGraph{N,E}() where {N<:AbstractNode,E<:UndirectedEdge}
        new(Dict(), Dict(), Dict())
    end
end

function MapUDGraph(nodes, edges)
    graph = MapUDGraph{Node,Edge}()
    for node in nodes
        updatenode!(graph, Node(), node)
    end
    for (i, edge) in enumerate(edges)
        updateedge!(graph, Edge(edge...), i)
    end
    graph
end


struct VectorUDGraph{N<:AbstractNode,E<:UndirectedEdge} <: UndirectedGraph
    nodes::Vector{N}
    edges::Vector{E}
    adjacency::Vector{Dict{Int,Int}}
    property::GraphPropertyVectors
end

function VectorUDGraph(size::Int, edges::AbstractArray{Tuple{Int,Int}})
    # do not use `fill`
    ns = [Node() for i in 1:size]
    adj = [Dict() for i in 1:size]
    es = []
    for (i, (u, v)) in enumerate(edges)
        push!(es, Edge(u, v))
        adj[u][v] = i
        adj[v][u] = i
    end
    VectorUDGraph{Node,Edge}(ns, es, adj, GraphPropertyVectors())
end

function VectorUDGraph{N,E}(graph::MapUDGraph{N,E}
        ) where {N<:AbstractNode,E<:UndirectedEdge}
    ns = []
    es = []
    adj = [Dict() for n in graph.nodes]
    nodemap = Dict()
    edgemap = Dict()
    # The order of node indices should be kept for some cannonicalization
    # operations (ex. chirality flag).
    nkeys = sort(collect(keys(graph.nodes)))
    ekeys = sort(collect(keys(graph.edges)))
    for (i, k) in enumerate(nkeys)
        nodemap[k] = i
        push!(ns, graph.nodes[k])
    end
    for (i, k) in enumerate(ekeys)
        edgemap[k] = i
        e = graph.edges[k]
        e.u = nodemap[e.u]
        e.v = nodemap[e.v]
        push!(es, e)
    end
    for (u, nbrs) in graph.adjacency
        for (v, e) in nbrs
            adj[nodemap[u]][nodemap[v]] = edgemap[e]
        end
    end
    VectorUDGraph{N,E}(ns, es, adj, GraphPropertyVectors())
end


getnode(graph::UndirectedGraph, idx) = graph.nodes[idx]

getedge(graph::UndirectedGraph, idx) = graph.edges[idx]
getedge(graph::UndirectedGraph, u, v) = getedge(graph, graph.adjacency[u][v])

hasedge(graph::UndirectedGraph, u, v) = haskey(graph.adjacency[u], v)

# TODO: `enumerate` yields `Tuple` whereas `Dict` yields `Pair`
nodesiter(graph::VectorUDGraph) = enumerate(graph.nodes)
nodesiter(graph::MapUDGraph) = graph.nodes
nodevector(graph::VectorUDGraph) = graph.nodes


nodekeys(graph::VectorUDGraph) = collect(1:nodecount(graph))
nodekeys(graph::MapUDGraph) = collect(keys(graph.nodes))

nodeset(graph::VectorUDGraph) = Set(1:nodecount(graph))
nodeset(graph::MapUDGraph) = Set(keys(graph.nodes))


# TODO: `enumerate` yields `Tuple` whereas `Dict` yields `Pair`
edgesiter(graph::VectorUDGraph) = enumerate(graph.edges)
edgesiter(graph::MapUDGraph) = graph.edges
edgevector(graph::VectorUDGraph) = graph.edges


edgekeys(graph::VectorUDGraph) = collect(1:edgecount(graph))
edgekeys(graph::MapUDGraph) = collect(keys(graph.edges))

edgeset(graph::VectorUDGraph) = Set(1:edgecount(graph))
edgeset(graph::MapUDGraph) = Set(keys(graph.edges))


neighbors(graph::UndirectedGraph, idx) = graph.adjacency[idx]


function updatenode!(graph::MapUDGraph, node, idx)
    """Add or update a node"""
    graph.nodes[idx] = node
    if !(idx in keys(graph.adjacency))
        graph.adjacency[idx] = Dict()
    end
    return
end


function updateedge!(graph::MapUDGraph, edge, idx)
    """Add or update an edge"""
    nodes = nodekeys(graph)
    (edge.u in nodes) || throw(KeyError(edge.u))
    (edge.v in nodes) || throw(KeyError(edge.v))
    graph.edges[idx] = edge
    graph.adjacency[edge.u][edge.v] = idx
    graph.adjacency[edge.v][edge.u] = idx
    return
end

updateedge!(
    G::MapUDGraph, edge, u, v) = updateedge!(G, edge, graph.adjacency[u][v])


function unlinknode!(graph::MapUDGraph, idx)
    """Remove a node and its connecting edges"""
    (idx in nodekeys(graph)) || throw(KeyError(idx))
    for (n, nbr) in neighbors(graph, idx)
        delete!(graph.edges, nbr)
        delete!(graph.adjacency[n], idx)
    end
    delete!(graph.nodes, idx)
    delete!(graph.adjacency, idx)
    return
end


function unlinkedge!(graph::MapUDGraph, u, v)
    """Remove an edge"""
    nodes = nodekeys(graph)
    (u in nodes) || throw(KeyError(u))
    (v in nodes) || throw(KeyError(v))
    delete!(graph.edges, graph.adjacency[u][v])
    delete!(graph.adjacency[u], v)
    delete!(graph.adjacency[v], u)
    return
end

function unlinkedge!(graph::MapUDGraph, idx)
    """Remove an edge"""
    (idx in keys(graph.edges)) || throw(KeyError(idx))
    e = getedge(graph, idx)
    delete!(graph.edges, idx)
    delete!(graph.adjacency[e.u], e.v)
    delete!(graph.adjacency[e.v], e.u)
    return
end


nodetype(graph::MapUDGraph) = valtype(graph.nodes)
nodetype(graph::VectorUDGraph) = eltype(graph.nodes)

edgetype(graph::MapUDGraph) = valtype(graph.edges)
edgetype(graph::VectorUDGraph) = eltype(graph.edges)

function similargraph(graph::UndirectedGraph)
    N = nodetype(graph)
    E = edgetype(graph)
    MapUDGraph{N,E}()
end


newgraph(graph::MapUDGraph) = graph # TODO:

"""
    newgraph(graph::VectorUDGraph)

Generate map graph from vector graph.
"""
function newgraph(graph::VectorUDGraph)
    newg = similargraph(graph)
    for i in nodekeys(graph)
        newg.nodes[i] = graph.nodes[i]
        newg.adjacency[i] = graph.adjacency[i]
    end
    for i in edgeskeys(graph)
        newg.edges[i] = graph.edges[i]
    end
    return newg
end
