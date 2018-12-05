#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    triangles


function triangles(graph::AbstractUGraph)
    triads = Set()
    for n in nodekeys(graph)
        nkeys = neighborkeys(graph, n)
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
