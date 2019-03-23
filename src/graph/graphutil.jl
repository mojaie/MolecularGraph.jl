#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    clone,
    neighbors,
    adjacencies, successors, predecessors,
    adjacentnodes, successornodes, predecessornodes,
    incidences, in_incidences, outincidences,
    incidentedges, inedges, outedges,
    hasproperty

clone(node::AbstractNode) = node

neighbors(g::DirectedGraph, i) = merge(outneighbors(g, i), inneighbors(g, i))

adjacencies(g::AbstractGraph, i) = Set(keys(neighbors(g, i)))
successors(g::DirectedGraph, i) = Set(keys(outneighbors(g, i)))
predecessors(g::DirectedGraph, i) = Set(keys(inneighbors(g, i)))

adjacentnodes(g::AbstractGraph, i) = getnode.((g,), adjacencies(g, i))
successornodes(g::DirectedGraph, i) = getnode.((g,), successors(g, i))
predecessornodes(g::DirectedGraph, i) = getnode.((g,), predecessors(g, i))

incidences(g::AbstractGraph, i) = Set(values(neighbors(g, i)))
outincidences(g::DirectedGraph, i) = Set(values(outneighbors(g, i)))
in_incidences(g::DirectedGraph, i) = Set(values(inneighbors(g, i)))

incidentedges(g::AbstractGraph, i) = getedge.((g,), incidences(g, i))
outedges(g::DirectedGraph, i) = getedge.((g,), outincidences(g, i))
inedges(g::DirectedGraph, i) = getedge.((g,), in_incidences(g, i))

neighborcount(g::AbstractGraph, i) = length(neighbors(g, i))
outdegree(g::DirectedGraph, i) = length(outneighbors(g, i))
indegree(g::DirectedGraph, i) = length(inneighbors(g, i))

nodecount(g::Union{DiGraph,Graph}) = length(g.nodes)
edgecount(g::Union{DiGraph,Graph}) = length(g.edges)

hasproperty(g::AbstractGraph, property::Symbol) = (
    isdefined(g, property) && !isempty(getproperty(g.property, property)))
