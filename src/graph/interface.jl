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
    AbstractNode,
    AbstractEdge,
    neighborkeys,
    neighbornodes,
    neighboredgekeys,
    neighboredges,
    neighborcount,
    degree


abstract type AbstractUGraph end

abstract type UGraph <: AbstractUGraph end
abstract type UGraphView <: AbstractUGraph end

abstract type MapUGraph <: UGraph end
abstract type VectorUGraph <: UGraph end

abstract type AbstractNode end
abstract type AbstractEdge end


neighborkeys(graph, idx) = collect(keys(neighbors(graph, idx)))
neighbornodes(graph, idx) = getnode.((graph,), neighborkeys(graph, idx))
neighboredgekeys(graph, idx) = collect(values(neighbors(graph, idx)))
neighboredges(graph, idx) = getedge.((graph,), neighboredgekeys(graph, idx))
neighborcount(graph, idx) = length(neighbors(graph, idx))
degree = neighborcount
