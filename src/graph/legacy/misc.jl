#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    twocoloring, triangles


"""
    twocoloring(graph::UndirectedGraph)

Do 2-coloring and return two sets of nodes.
"""
function twocoloring(graph::UndirectedGraph)
    nodes = nodeset(graph)
    c1 = Set()
    c2 = Set()
    while !isempty(nodes)
        root = pop!(nodes)
        queue = [root]
        push!(c1, root)
        while !isempty(queue)
            i = popfirst!(queue)
            for adj in adjacencies(graph, i)
                if !(adj in c1 || adj in c2)
                    if i in c1
                        push!(c2, adj)
                        push!(queue, adj)
                    else
                        push!(c1, adj)
                        push!(queue, adj)
                    end
                elseif i in c1 && adj in c1
                    return  # No 2-coloring
                elseif i in c2 && adj in c2
                    return  # No 2-coloring
                end
            end
        end
        setdiff!(nodes, c1, c2)
    end
    return c1, c2
end


"""
    triangles(graph::UndirectedGraph) -> Set{Tuple{Int,Int,Int}}

Fast computation of finding triangle node sets in the graph.
"""
function triangles(graph::UndirectedGraph)
    # deprecated: better implementation in Graphs.localclustering!
    triads = Set()
    for n in nodeset(graph)
        adjs = collect(adjacencies(graph, n))
        length(adjs) < 2 && continue
        for (e1, e2) in combinations(length(adjs))
            u, v = (adjs[e1], adjs[e2])
            if hasedge(graph, u, v)
                push!(triads, tuple(sort([n, u, v])...))
            end
        end
    end
    return triads
end