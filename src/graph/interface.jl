#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    AbstractUGraph,
    UGraph,
    UGraphView,
    MapUGraph,
    VectorUGraph,
    AbstractDGraph,
    DGraph,
    DGraphView,
    MapDGraph,
    AbstractNode,
    AbstractEdge,
    AbstractDirectedEdge,
    Node,
    VF2State


# Undirected graph

abstract type AbstractUGraph end
abstract type UGraph <: AbstractUGraph end
abstract type UGraphView <: AbstractUGraph end


# Directed graph

abstract type AbstractDGraph end
abstract type DGraph <: AbstractDGraph end
abstract type DGraphView <: AbstractDGraph end


# Components

abstract type AbstractNode end
abstract type AbstractEdge end
abstract type AbstractDirectedEdge end


# Node

struct Node <: AbstractNode
    attr::Dict
end

Node() = Node(Dict())


# Graph algorithm states

abstract type VF2State end
