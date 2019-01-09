#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    pathgraph,
    completegraph


"""
    pathgraph(length::Int) -> MapUDGraph{Node,Edge}

Generate path graph with the given length.
"""
function pathgraph(length::Int)
    return MapUDGraph(1:length, (i, i + 1) for i in 1:(length-1))
end


"""
    completegraph(type::MapUDGraph, length::Int)

Generate complete graph with the given size.
"""
function completegraph(size::Int)
    if size < 0
        throw(ValueError("Invalid value: $(size)"))
    elseif size == 0
        return MapUDGraph{Node,Edge}()
    elseif size == 1
        return MapUDGraph([1], [])
    end
    return MapUDGraph(1:size, ((u, v) for (u, v) in combinations(1:size)))
end
