#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    Edge, MapGraph, VectorGraph,
    setnodes, mapgraph, vectorgraph,
    nodevector, edgevector


struct Edge <: UndirectedEdge
    u::Int
    v::Int
end

setnodes(edge::Edge, u, v) = Edge(u, v)


struct MapGraph{N<:AbstractNode,E<:UndirectedEdge} <: Graph
    nodes::Dict{Int,N}
    edges::Dict{Int,E}
    neighbormap::Dict{Int,Dict{Int,Int}}

    function MapGraph{N,E}() where {N<:AbstractNode,E<:UndirectedEdge}
        new(Dict(), Dict(), Dict())
    end
end

"""
    mapgraph(::Type{N}, ::Type{E}
        ) where {N<:AbstractNode,E<:UndirectedEdge} -> MapGraph{N,E}()

Generate empty map graph that has nodes and edges with the given types.
"""
mapgraph(::Type{N}, ::Type{E}
    ) where {N<:AbstractNode,E<:UndirectedEdge} = MapGraph{N,E}()

"""
    mapgraph(nodes, edges) -> MapGraph{Node,Edge}

Generate map graph that has given nodes and edges represented by the list of
node indices in integer and the list of pairs of node indices, respectively.
"""
function mapgraph(nodes, edges)
    graph = MapGraph{Node,Edge}()
    for node in nodes
        graph.nodes[node] = Node()
        graph.neighbormap[node] = Dict()
    end
    for (i, edge) in enumerate(edges)
        (u, v) = edge
        graph.edges[i] = Edge(u, v)
        graph.neighbormap[u][v] = i
        graph.neighbormap[v][u] = i
    end
    return graph
end

"""
    mapgraph(graph::UndirectedGraph; clone=false) -> MapGraph

Convert the given graph into a new `MapGraph`. The node type and edge type are
inherited from the original graph. If the given graph is `MapGraph`, return
a copy of the graph.

If you really need mutable nodes and edges, and want the graph with deepcopied
elements, implement `clone` and `setnodes` as deepcopy methods for nodes and
edges, respectively.
"""
function mapgraph(graph::UndirectedGraph)
    newg = mapgraph(nodetype(graph), edgetype(graph))
    for (i, node) in nodesiter(graph)
        newg.nodes[i] = clone(node)
        newg.neighbormap[i] = Dict()
    end
    for (i, edge) in edgesiter(graph)
        newg.edges[i] = setnodes(edge, edge.u, edge.v)
        newg.neighbormap[edge.u][edge.v] = i
        newg.neighbormap[edge.v][edge.u] = i
    end
    return newg
end


struct VectorGraph{N<:AbstractNode,E<:UndirectedEdge} <: Graph
    nodes::Vector{N}
    edges::Vector{E}
    neighbormap::Vector{Dict{Int,Int}}
    property::GraphPropertyVectors

    function VectorGraph{N,E}() where {N<:AbstractNode,E<:UndirectedEdge}
        new([], [], [], GraphPropertyVectors())
    end
end

"""
    vectorgraph(::Type{N}, ::Type{E}
        ) where {N<:AbstractNode,E<:UndirectedEdge} -> VectorGraph{N,E}()

Generate empty vector graph that have nodes and edges with the given types.
"""
vectorgraph(::Type{N}, ::Type{E}
    ) where {N<:AbstractNode,E<:UndirectedEdge} = VectorGraph{N,E}()

"""
    vectorgraph(nodes, edges) -> VectorGraph{Node,Edge}

Generate vector graph that have given nodes and edges represented by the list of
node indices in integer and the list of pairs of node indices, respectively.
"""
function vectorgraph(size::Int, edges)
    graph = VectorGraph{Node,Edge}()
    append!(graph.nodes, [Node() for i in 1:size])
    append!(graph.neighbormap, [Dict() for i in 1:size])
    for (i, (u, v)) in enumerate(edges)
        push!(graph.edges, Edge(u, v))
        graph.neighbormap[u][v] = i
        graph.neighbormap[v][u] = i
    end
    return graph
end

"""
    vectorgraph(graph::UndirectedGraph; clone=false) -> VectorGraph

Convert the given graph into a new `VectorGraph`. The node type and edge type
are inherited from the original graph.

Node indices are sorted in ascending order and are re-indexed. This behavior is
intended for some cannonicalization operations (ex. chirality flag).

If you really need mutable nodes and edges, and want the graph with deepcopied
elements, implement `clone` and `setnodes` as deepcopy methods for nodes and
edges, respectively.
"""
function vectorgraph(graph::UndirectedGraph)
    newg = vectorgraph(nodetype(graph), edgetype(graph))
    nkeys = sort(nodekeys(graph))
    ekeys = sort(edgekeys(graph))
    nmap = Dict{Int,Int}()
    for (i, n) in enumerate(nkeys)
        nmap[n] = i
        node = getnode(graph, n)
        push!(newg.nodes, clone(node))
        push!(newg.neighbormap, Dict())
    end
    for (i, e) in enumerate(ekeys)
        edge = getedge(graph, e)
        u = nmap[edge.u]
        v = nmap[edge.v]
        push!(newg.edges, setnodes(edge, u, v))
        newg.neighbormap[u][v] = i
        newg.neighbormap[v][u] = i
    end
    return newg
end


getnode(graph::Graph, idx) = graph.nodes[idx]
getedge(graph::Graph, idx) = graph.edges[idx]
getedge(graph::Graph, u, v) = getedge(graph, graph.neighbormap[u][v])
hasedge(graph::Graph, u, v) = haskey(graph.neighbormap[u], v)

# TODO: `enumerate` yields `Tuple` whereas `Dict` yields `Pair`
nodesiter(graph::VectorGraph) = enumerate(graph.nodes)
nodesiter(graph::MapGraph) = graph.nodes
nodevector(graph::VectorGraph) = graph.nodes
edgesiter(graph::VectorGraph) = enumerate(graph.edges)
edgesiter(graph::MapGraph) = graph.edges
edgevector(graph::VectorGraph) = graph.edges

nodekeys(graph::VectorGraph) = collect(1:nodecount(graph))
nodekeys(graph::MapGraph) = collect(keys(graph.nodes))
nodeset(graph::VectorGraph) = Set(1:nodecount(graph))
nodeset(graph::MapGraph) = Set(keys(graph.nodes))
edgekeys(graph::VectorGraph) = collect(1:edgecount(graph))
edgekeys(graph::MapGraph) = collect(keys(graph.edges))
edgeset(graph::VectorGraph) = Set(1:edgecount(graph))
edgeset(graph::MapGraph) = Set(keys(graph.edges))

neighbors(graph::Graph, idx) = graph.neighbormap[idx]


function updatenode!(graph::MapGraph, node, idx)
    graph.nodes[idx] = node
    if !haskey(graph.neighbormap, idx)
        graph.neighbormap[idx] = Dict()
    end
    return
end

function updatenode!(graph::MapGraph, node)
    i = maximum(nodeset(graph)) + 1
    graph.nodes[i] = node
    graph.neighbormap[i] = Dict()
    return
end

function updatenode!(graph::VectorGraph, node, idx)
    idx > nodecount(graph) && throw(KeyError(idx))
    graph.nodes[idx] = node
    return
end

function updatenode!(graph::VectorGraph, node)
    push!(graph.nodes, node)
    push!(graph.neighbormap, Dict())
    return
end


function updateedge!(graph::MapGraph, edge, idx)
    nodes = nodeset(graph)
    (edge.u in nodes) || throw(KeyError(edge.u))
    (edge.v in nodes) || throw(KeyError(edge.v))
    if haskey(graph.edges, idx)
        old = getedge(graph, idx)
        delete!(graph.neighbormap[old.u], old.v)
        delete!(graph.neighbormap[old.v], old.u)
    end
    graph.edges[idx] = edge
    graph.neighbormap[edge.u][edge.v] = idx
    graph.neighbormap[edge.v][edge.u] = idx
    return
end

function updateedge!(graph::VectorGraph, edge, idx)
    idx > edgecount(graph) + 1 && throw(DomainError(idx))
    if idx == edgecount(graph) + 1
        updateedge!(graph, edge)
        return
    end
    ncnt = nodecount(graph)
    edge.u > ncnt && throw(DomainError(edge.u))
    edge.v > ncnt && throw(DomainError(edge.v))
    old = getedge(graph, idx)
    delete!(graph.neighbormap[old.u], old.v)
    delete!(graph.neighbormap[old.v], old.u)
    graph.edges[idx] = edge
    graph.neighbormap[edge.u][edge.v] = idx
    graph.neighbormap[edge.v][edge.u] = idx
    return
end

updateedge!(graph::MapGraph, edge
    ) = updateedge!(graph, edge, maximum(edgeset(graph)) + 1)

function updateedge!(graph::VectorGraph, edge)
    ncnt = nodecount(graph)
    edge.u > ncnt && throw(DomainError(edge.u))
    edge.v > ncnt && throw(DomainError(edge.v))
    i = edgecount(graph) + 1
    push!(graph.edges, edge)
    graph.neighbormap[edge.u][edge.v] = i
    graph.neighbormap[edge.v][edge.u] = i
    return
end

updateedge!(graph::Graph, edge, u, v
    ) = updateedge!(graph, edge, graph.neighbormap[u][v])



function unlinknode!(graph::MapGraph, idx)
    (idx in nodekeys(graph)) || throw(KeyError(idx))
    for (n, nbr) in neighbors(graph, idx)
        delete!(graph.edges, nbr)
        delete!(graph.neighbormap[n], idx)
    end
    delete!(graph.nodes, idx)
    delete!(graph.neighbormap, idx)
    return
end


function unlinkedge!(graph::MapGraph, u, v)
    nodes = nodekeys(graph)
    (u in nodes) || throw(KeyError(u))
    (v in nodes) || throw(KeyError(v))
    delete!(graph.edges, graph.neighbormap[u][v])
    delete!(graph.neighbormap[u], v)
    delete!(graph.neighbormap[v], u)
    return
end

function unlinkedge!(graph::MapGraph, idx)
    (idx in keys(graph.edges)) || throw(KeyError(idx))
    e = getedge(graph, idx)
    delete!(graph.edges, idx)
    delete!(graph.neighbormap[e.u], e.v)
    delete!(graph.neighbormap[e.v], e.u)
    return
end


nodetype(graph::MapGraph) = valtype(graph.nodes)
nodetype(graph::VectorGraph) = eltype(graph.nodes)

edgetype(graph::MapGraph) = valtype(graph.edges)
edgetype(graph::VectorGraph) = eltype(graph.edges)
