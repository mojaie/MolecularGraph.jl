#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    ancestors,
    descendants,
    topologicalsort


function ancestors(graph::DirectedGraph, idx::Int)
    rev = reversegraph(graph)
    return reachablenodes(rev, idx)
end


function descendants(graph::DirectedGraph, idx::Int)
    return reachablenodes(graph, idx)
end


function topologicalsort(graph::DirectedGraph)
    result = Int[]
    stack = [i for i in nodekeys(graph) if indegree(graph, i) == 0]
    used = Set{Int}()
    while !isempty(stack)
        n = pop!(stack)
        push!(result, n)
        for (s, edge) in outneighbors(graph, n)
            if edge in used
                continue
            end
            push!(used, edge)
            if isempty(setdiff(in_incidences(graph, s), used))
                push!(stack, s)
            end
        end
    end
    if edgecount(graph) != length(used)
        throw(ErrorException("cycle found"))
    end
    return result
end
