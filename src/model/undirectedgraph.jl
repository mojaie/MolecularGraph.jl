#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    Node,
    Edge,
    UndirectedGraph,
    getnode,
    getedge,
    neighbors,
    newnode!,
    updateedge!,
    unlinknode!,
    unlinkedge!


abstract type Node end
abstract type Edge end

mutable struct IntNode <: Node
    index::UInt32
end

mutable struct IntEdge <: Edge
    u::UInt32
    v::UInt32
end


mutable struct UndirectedGraph{T<:Unsigned}
    # TODO: use weakref map
    # TODO: parametric Index
    nodes::Vector{Node}
    edges::Vector{Edge}
    adjacency::Vector{Dict{T, Edge}}
    nodemap::Dict{T, Node}
    adjmap::Dict{T, Dict{T, Edge}}
end

UndirectedGraph{T}() where T <: Integer =
    UndirectedGraph{T}([], [], [], Dict(), Dict())

function UndirectedGraph(nodes::AbstractArray{T},
                         edges::AbstractArray{Tuple{T, T}}) where T <: Integer
    graph = UndirectedGraph{UInt32}()
    for node in nodes
        newnode!(graph, IntNode(node))
    end
    for edge in edges
        updateedge!(graph, IntEdge(edge...))
    end
    graph
end


function getnode(graph::UndirectedGraph, idx)
    graph.nodemap[idx]
end


function getedge(graph::UndirectedGraph, u, v)
    graph.adjmap[u][v]
end


function neighbors(graph::UndirectedGraph, idx)
    graph.adjmap[idx]
end


function newnode!(graph::UndirectedGraph{T}, node::Node) where T <: Unsigned
    """Add or init a node
    Adjacency of the node is initialized
    Old Node and adjacency object should be removed by clean!
    """
    push!(graph.nodes, node)
    graph.nodemap[node.index] = node
    adj = Dict{T, Edge}()
    push!(graph.adjacency, adj)
    graph.adjmap[node.index] = adj
    return
end


function updateedge!(graph::UndirectedGraph, edge::Edge)
    """Add or update an edge
    Old Edge object should be removed by clean!
    """
    push!(graph.edges, edge)
    graph.adjmap[edge.u][edge.v] = edge
    graph.adjmap[edge.v][edge.u] = edge
    return
end


function unlinknode!(graph::UndirectedGraph, idx)
    """Remove a node from mapping
    Old Node object should be removed by clean!
    """
    deleteat!(graph.nodemap, idx)
    deleteat!(graph.adjacency, idx)
    return
end


function unlinkedge!(graph::UndirectedGraph, u, v)
    """Remove an edge from mapping
    Old Edge object should be removed by clean!
    """
    deleteat!(graph.adjacency[u], v)
    deleteat!(graph.adjacency[v], u)
    return
end
