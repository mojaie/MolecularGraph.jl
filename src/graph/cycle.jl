#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    minimumcyclebasis



function findmincycle(graph, S)
    G = edgesubgraph(graph, setdiff(edgeset(graph), S))
    U = disjointunion(G, G)
    for s in S
        (u, v) = getedge(graph, s)
        u1 = getunionnode(U, 1, u)
        v1 = getunionnode(U, 1, v)
        u2 = getunionnode(U, 2, u)
        v2 = getunionnode(U, 2, v)
        addedge!(U, u1, v2, DisjointUnionEdge(1, s))
        addedge!(U, u2, v1, DisjointUnionEdge(1, s))
    end
    paths = []
    for n in nodeset(graph)
        n1 = getunionnode(U, 1, n)
        n2 = getunionnode(U, 2, n)
        sp = shortestpathedges(U, n1, n2)
        push!(paths, sp)
    end
    minpath = sortstablemin(paths, by=length)
    return [getsourceedge(U, e).sourcekey for e in minpath]
end


"""
    minimumcyclebasis(graph::UndirectedGraph) -> Vector{Vector{Int}}

Returns minimum cycle basis represented as an array of edge sequence that make up a cycle.

The node sequences can be obtained by using [`minimumcyclenodes`](@ref).

# Reference

1. de Pina, J.C. Applications of shortest path methods. PhD thesis, University of Amsterdam,
   Nether-lands (1995)
1. Kavitha, T., Mehlhorn, K., Michail, D. & Paluch, K. E. An ``\\tilde{O}(m^{2}n)`` Algorithm for
   Minimum Cycle Basis of Graphs. Algorithmica 52, 333–349 (2008).
"""
@cachefirst function minimumcyclebasis(graph::UndirectedGraph)
    cycles = Vector{Int}[]
    for conn in connectedcomponents(graph)
        subg = nodesubgraph(graph, conn)
        S = [Set(e) for e in cotree_edges(subg)]
        N = length(S)  # N: circuit rank
        for k in 1:N
            minedges = findmincycle(subg, S[k])
            push!(cycles, minedges)
            for i in (k + 1):N
                if length(intersect(S[i], minedges)) % 2 == 1
                    S[i] = symdiff(S[i], S[k])
                end
            end
        end
    end
    return cycles
end
