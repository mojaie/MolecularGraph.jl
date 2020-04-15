#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    triangles


function triangles(graph::UndirectedGraph)
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
