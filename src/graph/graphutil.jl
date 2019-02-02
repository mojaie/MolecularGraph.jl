#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    neighbors,
    neighborkeys, succkeys, predkeys,
    neighbornodes, succnodes, prednodes,
    neighboredgekeys, inedgekeys, outedgekeys,
    neighboredges, inedges, outedges,
    neighborcount, degree, indegree, outdegree,
    nodecount, edgecount


neighbors(g::DGraph, i) = merge(successors(g, i), predecessors(g, i))

neighborkeys(g::AbstractGraph, i) = Set(keys(neighbors(g, i)))
succkeys(g::DGraph, i) = Set(keys(successors(g, i)))
predkeys(g::DGraph, i) = Set(keys(predecessors(g, i)))

neighbornodes(g::AbstractGraph, i) = getnode.((g,), neighborkeys(g, i))
succnodes(g::DGraph, i) = getnode.((g,), succkeys(g, i))
prednodes(g::DGraph, i) = getnode.((g,), predkeys(g, i))

neighboredgekeys(g::AbstractGraph, i) = Set(values(neighbors(g, i)))
outedgekeys(g::DGraph, i) = Set(values(successors(g, i)))
inedgekeys(g::DGraph, i) = Set(values(predecessors(g, i)))

neighboredges(g::AbstractGraph, i) = getedge.((g,), neighboredgekeys(g, i))
outedges(g::DGraph, i) = getedge.((g,), outedgekeys(g, i))
inedges(g::DGraph, i) = getedge.((g,), inedgekeys(g, i))

neighborcount(g::AbstractGraph, i) = length(neighbors(g, i))
degree = neighborcount
outdegree(g::DGraph, i) = length(successors(g, i))
indegree(g::DGraph, i) = length(predecessors(g, i))

nodecount(g::Union{DirectedGraph,UndirectedGraph}) = length(g.nodes)
edgecount(g::Union{DirectedGraph,UndirectedGraph}) = length(g.edges)
