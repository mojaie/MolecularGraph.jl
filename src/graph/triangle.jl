#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    triangles


function triangles(graph::UDGraph)
    triads = Set()
    for n in nodekeys(graph)
        nkeys = neighborset(graph, n)
        if length(nkeys) < 2
            continue
        end
        subg = nodesubgraph(graph, nkeys)
        for (i, e) in edgesiter(subg)
            push!(triads, tuple(sort([n, e.u, e.v])...))
        end
    end
    return triads
end
