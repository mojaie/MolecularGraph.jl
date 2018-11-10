#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    LineGraphNode,
    LineGraphEdge,
    linegraph


struct LineGraphNode <: AbstractNode
    n1::Int
    n2::Int
end


struct LineGraphEdge <: AbstractEdge
    u::Int
    v::Int
    common::Int
end


function linegraph(G)
    L = MutableUDGraph{LineGraphNode, LineGraphEdge}()
    for (i, edge) in enumerate(G.edges)
        L.nodes[i] = LineGraphNode(edge.u, edge.v)
        L.adjacency[i] = Dict()
    end
    ecnt = 1
    for (i, adj) in enumerate(G.adjacency)
        if degree(G, i) <= 1
            continue
        end
        # TODO: is there a stuff like itertools.combination?
        for e1 in values(adj)
            for e2 in values(adj)
                if e1 < e2
                    L.edges[ecnt] = LineGraphEdge(e1, e2, i)
                    L.adjacency[e1][e2] = ecnt
                    L.adjacency[e2][e1] = ecnt
                    ecnt += 1
                end
            end
        end
    end
    L
end
