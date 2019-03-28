#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    nodevalues, edgevalues


GView = Union{GraphView,DiGraphView}

getnode(view::GView, idx) = getnode(view.graph, idx)
getedge(view::GView, idx) = getedge(view.graph, idx)
getedge(view::GView, u, v) = getedge(view.graph, u, v)
hasedge(view::GView, u, v) = hasedge(view.graph, u, v)

nodesiter(view::GView) = nodesiter(view.graph)
edgesiter(view::GView) = edgesiter(view.graph)
# TODO: for only VectorGraph
nodevalues(view::GView) = nodevalues(view.graph)
edgevalues(view::GView) = edgevalues(view.graph)

nodekeys(view::GView) = nodekeys(view.graph)
edgekeys(view::GView) = edgekeys(view.graph)
nodeset(view::GView) = nodeset(view.graph)
edgeset(view::GView) = edgeset(view.graph)

neighbors(view::GView, idx) = neighbors(view.graph, idx)
outneighbors(view::DiGraphView, idx) = outneighbors(view.graph, idx)
inneighbors(view::DiGraphView, idx) = inneighbors(view.graph, idx)

nodecount(view::GView) = nodecount(view.graph)
edgecount(view::GView) = edgecount(view.graph)

nodetype(view::GView) = nodetype(view.graph)
edgetype(view::GView) = edgetype(view.graph)

updatenode!(view::GView, node, i) = updatenode!(view.graph, node, i)
updatenode!(view::GView, node) = updatenode!(view.graph, node)
updateedge!(view::GView, edge, i) = updateedge!(view.graph, edge, i)
updateedge!(view::GView, edge) = updateedge!(view.graph, edge)
updateedge!(view::GView, edge, u, v) = updateedge!(view.graph, edge, u, v)
unlinknode!(view::GView, i) = unlinknode!(view.graph, i)
unlinkedge!(view::GView, e) = unlinkedge!(view.graph, e)
unlinkedge!(view::GView, u, v) = unlinkedge!(view.graph, u, v)
