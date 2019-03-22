#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    nodevector, edgevector,
    similargraph, newgraph


getnode(view::GraphView, idx) = getnode(view.graph, idx)

getedge(view::GraphView, idx) = getedge(view.graph, idx)
getedge(view::GraphView, u, v) = getedge(view.graph, u, v)

hasedge(view::GraphView, u, v) = hasedge(view.graph, u, v)

nodesiter(view::GraphView) = nodesiter(view.graph)
edgesiter(view::GraphView) = edgesiter(view.graph)
# TODO: for only VectorGraph
nodevector(view::GraphView) = nodevector(view.graph)
edgevector(view::GraphView) = edgevector(view.graph)

nodekeys(view::GraphView) = nodekeys(view.graph)
edgekeys(view::GraphView) = edgekeys(view.graph)

nodeset(view::GraphView) = nodeset(view.graph)
edgeset(view::GraphView) = edgeset(view.graph)

neighbors(view::GraphView, idx) = neighbors(view.graph, idx)

successors(view::DirectedGraphView, idx) = successors(view.graph, idx)
predecessors(view::DirectedGraphView, idx) = predecessors(view.graph, idx)

nodecount(view::GraphView) = nodecount(view.graph)
edgecount(view::GraphView) = edgecount(view.graph)

nodetype(view::GraphView) = nodetype(view.graph)
edgetype(view::GraphView) = edgetype(view.graph)
similargraph(view::GraphView) = similargraph(view.graph)


"""
    newgraph(view::UndirectedGraphView)

Generate graph from `GraphView`
"""
function newgraph(view::UndirectedGraphView; deepcopy=false)
    newg = similargraph(view.graph)
    for (i, node) in nodesiter(view)
        newg.nodes[i] = node
        newg.adjacency[i] = Dict{Int,Int}()
        for (nbr, e) in neighbors(view, i)
            if !(nbr in keys(newg.adjacency))
                edge = getedge(view, e) # TODO: deepcopy?
                edge.u = i
                edge.v = nbr
                newg.edges[e] = edge
            end
            newg.adjacency[i][nbr] = e
        end
    end
    return newg
end


# TODO: MutableGraphView
updatenode!(view::GraphView, node, i) = updatenode!(view.graph, node, i)
updateedge!(view::GraphView, edge, i) = updateedge!(view.graph, edge, i)
updateedge!(view::GraphView, edge, u, v) = updateedge!(view.graph, edge, u, v)
unlinknode!(view::GraphView, i) = unlinknode!(view.graph, i)
unlinkedge!(view::GraphView, e) = unlinkedge!(view.graph, e)
unlinkedge!(view::GraphView, u, v) = unlinkedge!(view.graph, u, v)
