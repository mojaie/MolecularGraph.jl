#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    getnode, getedge,
    nodesiter, edgesiter,
    nodekeys, edgekeys,
    neighbors, successors, predecessors,
    nodecount, edgecount,
    nodetype, edgetype, similarmap,
    updatenode!, updateedge!,
    unlinknode!, unlinkedge!


getnode(view::GraphView, idx) = getnode(view.graph, idx)

getedge(view::GraphView, idx) = getedge(view.graph, idx)
getedge(view::GraphView, u, v) = getedge(view.graph, u, v)

nodesiter(view::GraphView) = nodesiter(view.graph)
edgesiter(view::GraphView) = edgesiter(view.graph)

nodekeys(view::GraphView) = nodekeys(view.graph)
edgekeys(view::GraphView) = nodekeys(view.graph)

neighbors(view::GraphView, idx) = neighbors(view.graph, idx)

successors(view::DirectedGraphView, idx) = successors(view.graph, idx)
predecessors(view::DirectedGraphView, idx) = predecessors(view.graph, idx)

nodecount(view::GraphView) = nodecount(view.graph)
edgecount(view::GraphView) = edgecount(view.graph)

nodetype(view::GraphView) = nodetype(view.graph)
edgetype(view::GraphView) = edgetype(view.graph)
similarmap(view::GraphView) = similarmap(view.graph)

# TODO: MutableGraphView
updatenode!(view::GraphView, node, i) = updatenode!(view.graph, node, i)
updateedge!(view::GraphView, edge, i) = updateedge!(view.graph, edge, i)
updateedge!(view::GraphView, edge, u, v) = updateedge!(view.graph, edge, u, v)
unlinknode!(view::GraphView, i) = unlinknode!(view.graph, i)
unlinkedge!(view::GraphView, e) = unlinkedge!(view.graph, e)
unlinkedge!(view::GraphView, u, v) = unlinkedge!(view.graph, u, v)
