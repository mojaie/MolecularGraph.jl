#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    dfstree, dfstree_edges, cotree_edges, circuitrank


"""
    dfstree(adjfunc, graph, root; rootvalue = 0) -> Dict

Return tree node => predecessor dict of the connected component.
"""
function dfstree(adjfunc, graph, root; rootvalue = 0)
    stack = [root]
    pred = Dict(root => rootvalue)
    while !isempty(stack)
        i = pop!(stack)
        for adj in adjfunc(graph, i)
            if !haskey(pred, adj)
                pred[adj] = i
                push!(stack, adj)
            end
        end
    end
    return pred
end


"""
    edgedfstree(nbrfunc, graph, root; rootvalue = 0) -> Dict

Return tree node => edge to the predecessor dict of the connected component.
"""
function edgedfstree(nbrfunc, graph, root; rootvalue = 0)
    stack = [root]
    edges = Dict(root => rootvalue)
    while !isempty(stack)
        i = pop!(stack)
        for (inc, adj) in nbrfunc(graph, i)
            if !haskey(edges, adj)
                edges[adj] = inc
                push!(stack, adj)
            end
        end
    end
    return edges
end


"""
    cotree_edges(graph::UndirectedGraph, root) -> Set{Int}

Return a set of co-tree edges of the connected graph.
"""
function cotree_edges(graph::UndirectedGraph)
    root = pop!(nodeset(graph))
    treeedges = Set(values(edgedfstree(neighbors, graph, root)))
    return setdiff(edgeset(graph), treeedges)
end


"""
    circuitrank(graph::UndirectedGraph) -> Int

Return circuit rank (the number of cycles) of the undirected graph.
"""
function circuitrank(graph::UndirectedGraph)
    nodecount(graph) == 0 && return 0
    hascache(graph, :mincycles) && return length(graph.cache[:mincycles])
    conncount = length(connectedcomponents(graph))
    return edgecount(graph) - nodecount(graph) + conncount
end
