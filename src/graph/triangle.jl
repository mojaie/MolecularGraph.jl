#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    triangles


function triangles(graph::UndirectedGraph)
    triads = Set()
    for n in nodeset(graph)
        nkeys = adjacencies(graph, n)
        length(nkeys) < 2 && continue
        for (u, v) in combinations(nkeys)
            if hasedge(graph, u, v)
                push!(triads, tuple(sort([n, u, v])...))
            end
        end
    end
    return triads
end
