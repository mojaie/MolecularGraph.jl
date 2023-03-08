#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    line_graph


"""
    line_graph(G::SimpleGraph) -> SimpleGraph, Dict, Dict

Generate line graph, reverse mapping lg(v) -> g(e) and shared nodes mapping lg(e) -> g(v)
"""
function line_graph(g::SimpleGraph{T}) where T
    l = SimpleGraph{T}(ne(g))
    edge_rank = Dict(e => i for (i, e) in enumerate(edges(g)))
    revmap = Dict(i => e for (i, e) in enumerate(edges(g)))
    sharednode = Dict{Edge{T},T}()
    for i in vertices(g)
        degree(g, i) < 2 && continue
        nbrs = neighbors(g, i)
        for (m, n) in combinations(length(nbrs))
            e = undirectededge(T,
                edge_rank[undirectededge(T, i, nbrs[m])],
                edge_rank[undirectededge(T, i, nbrs[n])]
            )
            add_edge!(l, e)
            sharednode[e] = i
        end
    end
    return l, revmap, sharednode
end