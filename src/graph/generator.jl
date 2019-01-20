#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    pathgraph,
    cyclegraph,
    completegraph


"""
    pathgraph(length::Int) -> MapUDGraph{Node,Edge}

Generate path graph with the given size.
"""
function pathgraph(size::Int)
    if size < 2
        throw(ValueError("Invalid value: $(size)"))
    end
    return MapUDGraph(1:size, (i, i + 1) for i in 1:(size-1))
end


"""
    cyclegraph(length::Int) -> MapUDGraph{Node,Edge}

Generate cycle graph with the given size.
"""
function cyclegraph(size::Int)
    if size < 3
        throw(ValueError("Invalid value: $(size)"))
    end
    edges = [(i, i + 1) for i in 1:(size-1)]
    push!(edges, (size, 1))
    return MapUDGraph(1:size, edges)
end


"""
    completegraph(type::MapUDGraph, length::Int) -> MapUDGraph{Node,Edge}

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
