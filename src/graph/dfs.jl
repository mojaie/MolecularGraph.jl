#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    dfstree, dfstree_edges, cotree_edges, circuitrank


function dfstree(adjfunc, graph, root)
    stack = [root]
    pred = Dict{Int,Int}(root => 0)
    while !isempty(stack)
        i = pop!(stack)
        for adj in adjfunc(graph, i)
            if !haskey(pred, adj)
                pred[adj] = i
                push!(stack, adj)
            end
        end
    end
    delete!(pred, root)
    return pred
end


function dfstree_edges(nbrfunc, graph, root)
    stack = [root]
    edges = Dict{Int,Int}(root => 0)
    while !isempty(stack)
        i = pop!(stack)
        for (inc, adj) in nbrfunc(graph, i)
            if !haskey(edges, adj)
                edges[adj] = inc
                push!(stack, adj)
            end
        end
    end
    delete!(edges, root)
    return edges
end


"""
    cotree_edges(graph::UndirectedGraph, root) -> Set{Int}

Return a set of co-tree edges of a spanning tree.
"""
function cotree_edges(graph::UndirectedGraph, root)
    edges = dfstree_edges(neighbors, graph, root)
    return setdiff(edgeset(graph), values(edges))
end


"""
    circuitrank(graph::UndirectedGraph) -> Int

Return circuit rank (the number of cycles) of the graph.
"""
function circuitrank(graph::UndirectedGraph)
    nodecount(graph) == 0 && return 0
    if isdefined(graph, :cache) && haskey(graph.cache, :mincycles)
        return length(graph.cache[:mincycles])
    end
    treesize = 0
    for conn in connectedcomponents(graph)
        root = iterate(conn)[1]
        edges = dfstree_edges(neighbors, graph, root)
        treesize += length(edges)
    end
    return edgecount(graph) - treesize
end
