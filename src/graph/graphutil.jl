#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    neighbors,
    neighborset, succset, predset,
    neighbornodes, succnodes, prednodes,
    neighboredgeset, inedgeset, outedgeset,
    neighboredges, inedges, outedges,
    hasproperty


neighbors(g::DGraph, i) = merge(successors(g, i), predecessors(g, i))

neighborset(g::AbstractGraph, i) = Set(keys(neighbors(g, i)))
succset(g::DGraph, i) = Set(keys(successors(g, i)))
predset(g::DGraph, i) = Set(keys(predecessors(g, i)))

neighbornodes(g::AbstractGraph, i) = getnode.((g,), neighborset(g, i))
succnodes(g::DGraph, i) = getnode.((g,), succkeys(g, i))
prednodes(g::DGraph, i) = getnode.((g,), predkeys(g, i))

neighboredgeset(g::AbstractGraph, i) = Set(values(neighbors(g, i)))
outedgeset(g::DGraph, i) = Set(values(successors(g, i)))
inedgeset(g::DGraph, i) = Set(values(predecessors(g, i)))

neighboredges(g::AbstractGraph, i) = getedge.((g,), neighboredgeset(g, i))
outedges(g::DGraph, i) = getedge.((g,), outedgeset(g, i))
inedges(g::DGraph, i) = getedge.((g,), inedgeset(g, i))

neighborcount(g::AbstractGraph, i) = length(neighbors(g, i))
outdegree(g::DGraph, i) = length(successors(g, i))
indegree(g::DGraph, i) = length(predecessors(g, i))

nodecount(g::Union{DirectedGraph,UndirectedGraph}) = length(g.nodes)
edgecount(g::Union{DirectedGraph,UndirectedGraph}) = length(g.edges)

hasproperty(g::AbstractGraph, property::Symbol) = (
    isdefined(g, property) && !isempty(getproperty(g.property, property)))
