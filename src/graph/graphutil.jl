#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    succkeys, predkeys, neighborkeys,
    succnodes, prednodes, neighbornodes,
    inedgekeys, outedgekeys, neighboredgekeys,
    inedges, outedges, neighboredges,
    indegree, outdegree, neighborcount, degree,
    nodecount, edgecount


succkeys(g::DGraph, i) = collect(keys(successors(g, i)))
predkeys(g::DGraph, i) = collect(keys(predecessors(g, i)))
neighborkeys(g::AbstractGraph, i) = collect(keys(neighbors(g, i)))

succnodes(g::DGraph, i) = getnode.((g,), succkeys(g, i))
prednodes(g::DGraph, i) = getnode.((g,), predkeys(g, i))
neighbornodes(g::AbstractGraph, i) = getnode.((g,), neighborkeys(g, i))

outedgekeys(g::DGraph, i) = collect(values(successors(g, i)))
inedgekeys(g::DGraph, i) = collect(values(predecessors(g, i)))
neighboredgekeys(g::AbstractGraph, i) = collect(values(neighbors(g, i)))

outedges(g::DGraph, i) = getedge.((g,), outedgekeys(g, i))
inedges(g::DGraph, i) = getedge.((g,), inedgekeys(g, i))
neighboredges(g::AbstractGraph, i) = getedge.((g,), neighboredgekeys(g, i))

outdegree(g::DGraph, i) = length(successors(g, i))
indegree(g::DGraph, i) = length(predecessors(g, i))
neighborcount(g::AbstractGraph, i) = length(neighbors(g, i))
degree = neighborcount

nodecount(g::Union{DirectedGraph,UndirectedGraph}) = length(g.nodes)
edgecount(g::Union{DirectedGraph,UndirectedGraph}) = length(g.edges)
