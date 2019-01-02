#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    ancestors,
    descendants,
    topologicalsort


function ancestors(graph::DGraph, idx::Int)
    rev = ReverseGraph(graph)
    return reachablenodes(rev, idx)
end


function descendants(graph::DGraph, idx::Int)
    return reachablenodes(graph, idx)
end


function topologicalsort(graph::DGraph)
    result = Int[]
    stack = [i for i in nodekeys(graph) if indegree(graph, i) == 0]
    used = Set{Int}()
    while !isempty(stack)
        n = pop!(stack)
        push!(result, n)
        for (s, edge) in successors(graph, n)
            if edge in used
                continue
            end
            push!(used, edge)
            if isempty(setdiff(inedgekeys(graph, s), used))
                push!(stack, s)
            end
        end
    end
    if edgecount(graph) != length(used)
        throw(OperationError("Cycle found"))
    end
    return result
end
