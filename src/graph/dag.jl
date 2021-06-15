#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    descendants, ancestors,
    topologicalsort


descendants(graph::DirectedGraph, idx::Int) = reachablenodes(graph, idx)
ancestors(graph::DirectedGraph, idx::Int) = reversereachablenodes(graph, idx)


"""
    topologicalsort(graph::DirectedGraph) -> Vector{Int}

Run topological sort to obtain an array of sorted nodes.
"""
function topologicalsort(graph::DirectedGraph)
    result = Int[]
    stack = [i for i in nodeset(graph) if indegree(graph, i) == 0]
    used = Set{Int}()
    while !isempty(stack)
        n = pop!(stack)
        push!(result, n)
        for (ninc, nadj) in outneighbors(graph, n)
            ninc in used && continue
            push!(used, ninc)
            if isempty(setdiff(in_incidences(graph, nadj), used))
                push!(stack, nadj)
            end
        end
    end
    edgecount(graph) == length(used) || throw(ErrorException("cycle found"))
    return result
end
