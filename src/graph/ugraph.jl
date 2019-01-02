#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    Edge,
    MapUGraph,
    VectorUGraph,
    connect,
    getnode,
    getedge,
    nodesiter,
    edgesiter,
    nodekeys,
    edgekeys,
    neighbors,
    nodecount,
    edgecount,
    neighborkeys,
    neighbornodes,
    neighboredgekeys,
    neighboredges,
    neighborcount,
    degree,
    updatenode!,
    updateedge!,
    unlinknode!,
    unlinkedge!,
    similarmap


struct Edge <: AbstractEdge
    u::Int
    v::Int
    attr::Dict
end

Edge(u, v) = Edge(u, v, Dict())
connect(e::Edge, u, v) = Edge(u, v, e.attr)


struct MapUGraph{N<:AbstractNode,E<:AbstractEdge} <: UGraph
    nodes::Dict{Int,N}
    edges::Dict{Int,E}
    adjacency::Dict{Int,Dict{Int,Int}}

    function MapUGraph{N,E}() where {N<:AbstractNode,E<:AbstractEdge}
        new(Dict(), Dict(), Dict())
    end
end

function MapUGraph(nodes::AbstractArray{Int},
                        edges::AbstractArray{Tuple{Int,Int}})
    graph = MapUGraph{Node,Edge}()
    for node in nodes
        updatenode!(graph, Node(), node)
    end
    for (i, edge) in enumerate(edges)
        updateedge!(graph, Edge(edge...), i)
    end
    graph
end


struct VectorUGraph{N<:AbstractNode,E<:AbstractEdge} <: UGraph
    nodes::Vector{N}
    edges::Vector{E}
    adjacency::Vector{Dict{Int,Int}}
end

function VectorUGraph(size::Int, edges::AbstractArray{Tuple{Int,Int}})
    # do not use `fill`
    ns = [Node() for i in 1:size]
    adj = [Dict() for i in 1:size]
    es = []
    for (i, (u, v)) in enumerate(edges)
        push!(es, Edge(u, v))
        adj[u][v] = i
        adj[v][u] = i
    end
    VectorUGraph{Node,Edge}(ns, es, adj)
end

function VectorUGraph{N,E}(graph::MapUGraph{N,E}
        ) where {N<:AbstractNode,E<:AbstractEdge}
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
        push!(es, connect(e, nodemap[e.u], nodemap[e.v]))
    end
    for (u, nbrs) in graph.adjacency
        for (v, e) in nbrs
            adj[nodemap[u]][nodemap[v]] = edgemap[e]
        end
    end
    VectorUGraph{N,E}(ns, es, adj)
end


getnode(graph::UGraph, idx) = graph.nodes[idx]

getedge(graph::UGraph, idx) = graph.edges[idx]
getedge(graph::UGraph, u, v) = getedge(graph, graph.adjacency[u][v])

nodesiter(graph::VectorUGraph) = enumerate(graph.nodes)
nodesiter(graph::MapUGraph) = graph.nodes

nodekeys(graph::VectorUGraph) = Set(1:nodecount(graph))
nodekeys(graph::MapUGraph) = Set(keys(graph.nodes))

edgesiter(graph::VectorUGraph) = enumerate(graph.edges)
edgesiter(graph::MapUGraph) = graph.edges

edgekeys(graph::VectorUGraph) = Set(1:edgecount(graph))
edgekeys(graph::MapUGraph) = Set(keys(graph.edges))

neighbors(graph::UGraph, idx) = graph.adjacency[idx]

nodecount(graph::UGraph) = length(graph.nodes)
edgecount(graph::UGraph) = length(graph.edges)

neighborkeys(graph::AbstractUGraph, idx) = collect(keys(neighbors(graph, idx)))
neighbornodes(
    graph::AbstractUGraph, idx) = getnode.((graph,), neighborkeys(graph, idx))
neighboredgekeys(
    graph::AbstractUGraph, idx) = collect(values(neighbors(graph, idx)))
neighboredges(
    graph::AbstractUGraph, idx
) = getedge.((graph,), neighboredgekeys(graph, idx))
neighborcount(
    graph::AbstractUGraph, idx) = length(neighbors(graph, idx))
degree = neighborcount


function updatenode!(graph::MapUGraph, node, idx)
    """Add or update a node"""
    graph.nodes[idx] = node
    if !(idx in keys(graph.adjacency))
        graph.adjacency[idx] = Dict()
    end
    return
end


function updateedge!(graph::MapUGraph, edge, idx)
    """Add or update an edge"""
    if !(edge.u in keys(graph.nodes))
        throw(OperationError("Missing node: $(edge.u)"))
    elseif !(edge.v in keys(graph.nodes))
        throw(OperationError("Missing node: $(edge.v)"))
    end
    graph.edges[idx] = edge
    graph.adjacency[edge.u][edge.v] = idx
    graph.adjacency[edge.v][edge.u] = idx
    return
end

updateedge!(
    G::MapUGraph, edge, u, v) = updateedge!(G, edge, graph.adjacency[u][v])


function unlinknode!(graph::MapUGraph, idx)
    """Remove a node and its connecting edges"""
    if !(idx in keys(graph.nodes))
        throw(OperationError("Missing node: $(idx)"))
    end
    for (n, nbr) in graph.adjacency[idx]
        delete!(graph.edges, nbr)
        delete!(graph.adjacency[n], idx)
    end
    delete!(graph.nodes, idx)
    delete!(graph.adjacency, idx)
    return
end


function unlinkedge!(graph::MapUGraph, u, v)
    """Remove an edge"""
    if !(u in keys(graph.nodes))
        throw(OperationError("Missing node: $(u)"))
    elseif !(v in keys(graph.nodes))
        throw(OperationError("Missing node: $(v)"))
    end
    delete!(graph.edges, graph.adjacency[u][v])
    delete!(graph.adjacency[u], v)
    delete!(graph.adjacency[v], u)
    return
end

function unlinkedge!(graph::MapUGraph, idx)
    """Remove an edge"""
    if !(idx in keys(graph.edges))
        throw(OperationError("Missing edge: $(idx)"))
    end
    e = getedge(graph, idx)
    delete!(graph.edges, idx)
    delete!(graph.adjacency[e.u], e.v)
    delete!(graph.adjacency[e.v], e.u)
    return
end


function similarmap(graph::MapUGraph)
    N = valtype(graph.nodes)
    E = valtype(graph.edges)
    MapUGraph{N,E}()
end

function similarmap(graph::VectorUGraph)
    N = eltype(graph.nodes)
    E = eltype(graph.edges)
    MapUGraph{N,E}()
end
