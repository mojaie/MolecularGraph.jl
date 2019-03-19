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
    nodekeys, nodeset, edgekeys, edgeset


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
    nodekeys(graph::AbstractGraph) -> Vector{Int}

Return graph node keys. If the given graph is a vector graph, the keys are in
index order, whereas the order of indices in map graph is not guaranteed.
"""
nodekeys


"""
    nodeset(graph::AbstractGraph) -> Set{Int}

Return the set of node keys.
"""
nodeset


"""
    edgekeys(graph::AbstractGraph) -> Vector{Int}

Return graph edge keys. If the given graph is a vector graph, the keys are in
index order, whereas the order of indices in map graph is not guaranteed.
"""
edgekeys


"""
    edgeset(graph::AbstractGraph) -> Set{Int}

Return the set of edge keys.
"""
edgeset
