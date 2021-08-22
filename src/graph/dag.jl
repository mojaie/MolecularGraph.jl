#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    descendants, ancestors, roots, leafs,
    topologicalsort, reversetopologicalsort,
    transitive_reduction


descendants(graph::DirectedGraph, idx::Int) = reachablenodes(graph, idx)
ancestors(graph::DirectedGraph, idx::Int) = reversereachablenodes(graph, idx)

roots(graph::OrderedDiGraph) = [i for i in 1:nodecount(graph) if indegree(graph, i) == 0]
roots(graph::DirectedGraph) = [i for i in nodeset(graph) if indegree(graph, i) == 0]

leafs(graph::OrderedDiGraph) = [i for i in 1:nodecount(graph) if outdegree(graph, i) == 0]
leafs(graph::DirectedGraph) = [i for i in nodeset(graph) if outdegree(graph, i) == 0]


"""
    topologicalsort(graph::DirectedGraph) -> Vector{Int}
    topologicalsort(graph::DirectedGraph, root) -> Vector{Int}

Run topological sort to obtain an array of sorted nodes.

If `root` is given, the method is applied to only the downstream nodes of `root`.
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



function topologicalsort(graph::DirectedGraph, root)
    sub = nodesubgraph(graph, union(descendants(graph, root), root))
    return topologicalsort(sub)
end


function reversetopologicalsort(graph::DirectedGraph)
    result = Int[]
    stack = [i for i in nodeset(graph) if outdegree(graph, i) == 0]
    used = Set{Int}()
    while !isempty(stack)
        n = pop!(stack)
        push!(result, n)
        for (ninc, nadj) in inneighbors(graph, n)
            ninc in used && continue
            push!(used, ninc)
            if isempty(setdiff(outincidences(graph, nadj), used))
                push!(stack, nadj)
            end
        end
    end
    edgecount(graph) == length(used) || throw(ErrorException("cycle found"))
    return result
end

function reversetopologicalsort(graph::DirectedGraph, leaf)
    sub = nodesubgraph(graph, union(ancestors(graph, leaf), leaf))
    return topologicalsort(sub)
end


"""
    longestpathedges(graph::DirectedGraph, root) -> ValueIterator

Return edges belong to the longest path from the given root node.
"""
function longestpathedges(graph::DirectedGraph, root)
    sub = nodesubgraph(graph, union(descendants(graph, root), root))
    dist = Dict{Int,Int}()
    inedges = Dict{Int,Int}()
    for n in nodeset(sub)
        dist[n] = n == root ? 0 : typemin(Int)
    end
    for n in topologicalsort(sub)
        for (inc, adj) in outneighbors(sub, n)
            if dist[adj] < dist[n] + 1
                dist[adj] = dist[n] + 1
                inedges[adj] = inc
            end
        end
    end
    return values(inedges)
end


"""
    transitive_reduction(graph::DirectedGraph) -> Set{Int}

Return transitive reduction as a subset of the DAG edges.
"""
function transitive_reduction(graph::DirectedGraph)
    edges = Set{Int}()
    for v in nodeset(graph)
        union!(edges, longestpathedges(graph, v))
    end
    return edges
end