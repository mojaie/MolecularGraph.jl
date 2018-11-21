#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    AbstractNode,
    AbstractEdge,
    AbstractUDGraph,
    Node,
    Edge,
    UDGraph,
    MutableUDGraph,
    connect,
    getnode,
    getedge,
    nodekeys,
    neighbors,
    neighborkeys,
    neighbornodes,
    neighboredgekeys,
    neighboredges,
    nodecount,
    edgecount,
    neighborcount,
    degree,
    updatenode!,
    updateedge!,
    unlinknode!,
    unlinkedge!


abstract type AbstractNode end
abstract type AbstractEdge end
abstract type AbstractUDGraph end


struct Node <: AbstractNode
    attr::Dict
end

Node() = Node(Dict())


struct Edge <: AbstractEdge
    u::Int
    v::Int
    attr::Dict
end

Edge(u, v) = Edge(u, v, Dict())
connect(e::Edge, u, v) = Edge(u, v, e.attr)


struct MutableUDGraph{N<:AbstractNode,E<:AbstractEdge} <: AbstractUDGraph
    nodes::Dict{Int,N}
    edges::Dict{Int,E}
    adjacency::Dict{Int,Dict{Int,Int}}

    function MutableUDGraph{N,E}() where {N<:AbstractNode,E<:AbstractEdge}
        new(Dict(), Dict(), Dict())
    end
end

function MutableUDGraph(nodes::AbstractArray{Int},
                        edges::AbstractArray{Tuple{Int,Int}})
    graph = MutableUDGraph{Node,Edge}()
    for node in nodes
        updatenode!(graph, Node(), node)
    end
    for (i, edge) in enumerate(edges)
        updateedge!(graph, Edge(edge...), i)
    end
    graph
end


struct UDGraph{N<:AbstractNode,E<:AbstractEdge} <: AbstractUDGraph
    nodes::Vector{N}
    edges::Vector{E}
    adjacency::Vector{Dict{Int,Int}}
end

function UDGraph(size::Int, edges::AbstractArray{Tuple{Int,Int}})
    # do not use `fill`
    ns = [Node() for i in 1:size]
    adj = [Dict() for i in 1:size]
    es = []
    for (i, (u, v)) in enumerate(edges)
        push!(es, Edge(u, v))
        adj[u][v] = i
        adj[v][u] = i
    end
    UDGraph{Node,Edge}(ns, es, adj)
end

function UDGraph{N,E}(graph::MutableUDGraph{N,E}
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
    UDGraph{N,E}(ns, es, adj)
end


getnode(graph::AbstractUDGraph, idx) = graph.nodes[idx]

getedge(graph::AbstractUDGraph, idx) = graph.edges[idx]
getedge(graph::AbstractUDGraph, u, v) = getedge(graph, graph.adjacency[u][v])


nodekeys(graph::UDGraph) = Set(1:nodecount(graph))
nodekeys(graph::MutableUDGraph) = Set(keys(graph.nodes))


neighbors(graph::AbstractUDGraph, idx) = graph.adjacency[idx]

neighborkeys(graph, idx) = collect(keys(graph.adjacency[idx]))

neighbornodes(graph, idx) = getnode.((graph,), keys(graph.adjacency[idx]))

neighboredgekeys(graph, idx) = collect(values(graph.adjacency[idx]))

neighboredges(graph, idx) = getedge.((graph,), values(graph.adjacency[idx]))


nodecount(graph::AbstractUDGraph) = length(graph.nodes)

edgecount(graph::AbstractUDGraph) = length(graph.edges)

neighborcount(graph, idx) = length(graph.adjacency[idx])
degree = neighborcount

function updatenode!(graph::MutableUDGraph, node, idx)
    """Add or update a node"""
    graph.nodes[idx] = node
    if !(idx in keys(graph.adjacency))
        graph.adjacency[idx] = Dict()
    end
    return
end


function updateedge!(graph::MutableUDGraph, edge, idx)
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

updateedge!(G, edge, u, v) = updateedge!(G, edge, graph.adjacency[u][v])


function unlinknode!(graph::MutableUDGraph, idx)
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


function unlinkedge!(graph::MutableUDGraph, u, v)
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

function unlinkedge!(graph::MutableUDGraph, idx)
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
