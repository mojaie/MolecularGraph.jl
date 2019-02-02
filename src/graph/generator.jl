#
# This file is a part of MolecularGraph.jl
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
        throw(DomainError(size, "path graph size should be 2 or more"))
    end
    return MapUDGraph(1:size, (i, i + 1) for i in 1:(size-1))
end


"""
    cyclegraph(length::Int) -> MapUDGraph{Node,Edge}

Generate cycle graph with the given size.
"""
function cyclegraph(size::Int)
    if size < 3
        throw(DomainError(size, "cycle graph size should be 3 or more"))
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
        throw(DomainError(size, "graph size should not be negative"))
    elseif size == 0
        return MapUDGraph{Node,Edge}()
    elseif size == 1
        return MapUDGraph([1], [])
    end
    return MapUDGraph(1:size, ((u, v) for (u, v) in combinations(1:size)))
end
