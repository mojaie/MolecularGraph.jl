#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    AbstractGraph,
    UndirectedGraph, DirectedGraph,
    UndirectedGraphView, DirectedGraphView,
    UDGraph, DGraph, GraphView,
    AbstractNode,
    UndirectedEdge, DirectedEdge,
    Node, GraphPropertyVectors,
    getnode, getedge, hasedge,
    nodesiter, egdesiter, nodekeys, edgekeys, nodeset, edgeset,
    neighbors, successors, predecessors,
    updatenode!, updateedge!, unlinknode!, unlinkedge!,
    nodecount, edgecount, neighborcount, degree, indegree, outdegree,
    nodetype, edgetype


abstract type AbstractGraph end
abstract type UndirectedGraph <: AbstractGraph end
abstract type DirectedGraph <: AbstractGraph end
abstract type UndirectedGraphView <: AbstractGraph end
abstract type DirectedGraphView <: AbstractGraph end

# Union types
# TODO: use traits
# https://github.com/JuliaLang/julia/issues/2345

UDGraph = Union{UndirectedGraph,UndirectedGraphView}
DGraph = Union{DirectedGraph,DirectedGraphView}
GraphView = Union{DirectedGraphView,UndirectedGraphView}


# Components

abstract type AbstractNode end
abstract type AbstractEdge end
abstract type UndirectedEdge <: AbstractEdge end
abstract type DirectedEdge <: AbstractEdge end


# Node

mutable struct Node <: AbstractNode
    attr::Dict
end

Node() = Node(Dict())


# Precalculated properties

struct GraphPropertyVectors
    mincycles::Vector{Vector{Int}}
    nodescycles::Vector{Set{Int}}
    nodescyclesizes::Vector{Set{Int}}
    edgescycles::Vector{Set{Int}}
    edgescyclesizes::Vector{Set{Int}}

    connected::Vector{Vector{Int}}
    connmembership::Vector{Int}
    biconnected::Vector{Vector{Int}}
    biconnmembership::Vector{Union{Int,Nothing}}
    cutvertices::Vector{Int}
    bridges::Vector{Int}
    twoedge::Vector{Vector{Int}}
    twoedgemembership::Vector{Union{Int,Nothing}}

    GraphPropertyVectors() = new(
        [], [], [], [], [], [], [], [], [], [],
        [], [], [])
end


"""
    getnode(graph::AbstractGraph, index) -> AbstractNode

Retrieve the node object at the given index within the graph.
"""
function getnode end


"""
    getedge(graph::AbstractGraph, index) -> AbstractEdge

Retrieve the edge object at the given index within the graph.
"""
function getedge(graph::AbstractGraph, index) end


"""
    getedge(graph::AbstractGraph, u, v) -> AbstractEdge

Retrieve an edge object which connects the given nodes.
"""
function getedge(graph::AbstractGraph, u, v) end


"""
    hasedge(graph::AbstractGraph, u, v) -> AbstractEdge

Return whether the given two nodes are connected by at least one edge.
"""
function hasedge end


"""
    nodesiter(graph::AbstractGraph)

An iterator that yields `(i, n)` where `i` is the node index, and `n` is the
node object at the index `i` within the graph.
"""
function nodesiter end


"""
    edgesiter(graph::AbstractGraph)

An iterator that yields `(i, e)` where `i` is the edge index, and `e` is the
edge object at the index `i` within the graph.
"""
function edgesiter end


"""
    nodekeys(graph::AbstractGraph) -> Vector{Int}

Return graph node keys. If the given graph is a vector graph, the keys are in
ascending order, whereas the order of indices in map graph is not guaranteed.
"""
function nodekeys end


"""
    nodeset(graph::AbstractGraph) -> Set{Int}

Return the set of node keys.
"""
function nodeset end


"""
    edgekeys(graph::AbstractGraph) -> Vector{Int}

Return graph edge keys. If the given graph is a vector graph, the keys are in
ascending order, whereas the order of indices in map graph is not guaranteed.
"""
function edgekeys end


"""
    edgeset(graph::AbstractGraph) -> Set{Int}

Return the set of edge keys.
"""
function edgeset end


"""
    neighbors(graph::AbstractGraph, n) -> Dict{Int,Int}

Return the mapping of adjacent node keys and incident edge keys connected to
the given node. If the graph is directed graph, both successors and
predecessors are mapped.
"""
function neighbors end


"""
    successors(graph::DirectedGraph, n) -> Dict{Int,Int}

Return the mapping of successor node keys and out edge keys connected to
the given node.
"""
function successors end


"""
    predecessors(graph::DirectedGraph, n) -> Dict{Int,Int}

Return the mapping of predecessor node keys and in edge keys connected to
the given node.
"""
function predecessors end


"""
    updatenode!(graph::AbstractGraph, node::AbstractNode, n)

Rebind the node object stored at `n` of the graph by the given node object.
If the index does not exist, add new node to the position `n`.
"""
function updatenode! end


"""
    updateedge!(graph::AbstractGraph, edge::AbstractEdge, e)

Rebind the edge object stored at `e` of the graph by the given edge object.
If the index does not exist, add new edge to the position `e`.
"""
function updateedge!(graph::AbstractGraph, edge, e) end


"""
    updateedge!(graph::AbstractGraph, edge::AbstractEdge, u, v)

Rebind the edge that connects nodes `u` and `v` by the given edge object.
If the nodes do not exist, throws `KeyError`.
"""
function updateedge!(graph::AbstractGraph, edge, u, v) end


"""
    unlinknode!(graph::AbstractGraph, n)

Delete the node at the index of `n` and its incident edges.
"""
function unlinknode! end


"""
    unlinkedge!(graph::AbstractGraph, e)

Delete the edge at the index of `e`.
"""
function unlinkedge!(graph::AbstractGraph, e) end


"""
    unlinkedge!(graph::AbstractGraph, u, v)

Delete the edge that connect nodes `u` and `v`.
"""
function unlinkedge!(graph::AbstractGraph, u, v) end


"""
    nodecount(graph::AbstractGraph) -> Int

Return the number of graph nodes.
"""
function nodecount end


"""
    edgecount(graph::AbstractGraph) -> Int

Return the number of graph edges.
"""
function edgecount end


"""
    neighborcount(graph::AbstractGraph, n) -> Int
    degree(graph::AbstractGraph, n) -> Int

Return the number of adjacent nodes of the node 'n'.
"""
function neighborcount end
degree = neighborcount


"""
    indegree(graph::DirectedGraph, n) -> Int

Return the number of predecessors of the node 'n'.
"""
function indegree end


"""
    outdegree(graph::DirectedGraph, n) -> Int

Return the number of successors of the node 'n'.
"""
function outdegree end


"""
    nodetype(graph::AbstractGraph)

Return the node type of the graph
"""
function nodetype end

"""
    edgetype(graph::AbstractGraph)

Return the edge type of the graph
"""
function edgetype end
