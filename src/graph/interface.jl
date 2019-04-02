#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    AbstractGraph,
    UndirectedGraph, DirectedGraph,
    GraphView, DiGraphView,
    AbstractNode,
    UndirectedEdge, DirectedEdge,
    Node, GraphPropertyVectors,
    getnode, getedge, hasedge,
    nodekeys, edgekeys, nodevalues, edgevalues,
    nodesiter, edgesiter, nodeset, edgeset,
    neighbors, outneighbors, inneighbors,
    updatenode!, updateedge!, unlinknodes, unlinkedges,
    nodecount, edgecount, neighborcount, degree, indegree, outdegree,
    nodetype, edgetype,
    @cache, clearcache!


abstract type AbstractGraph end
abstract type UndirectedGraph <: AbstractGraph end
abstract type DirectedGraph <: AbstractGraph end
abstract type GraphView <: UndirectedGraph end
abstract type DiGraphView <: DirectedGraph end

# Union types
# TODO: use traits
# https://github.com/JuliaLang/julia/issues/2345


# Components

abstract type AbstractNode end
abstract type AbstractEdge end
abstract type UndirectedEdge <: AbstractEdge end
abstract type DirectedEdge <: AbstractEdge end


struct Node <: AbstractNode
    dummy::Nothing
end

Node() = Node(nothing)


"""
    getnode(graph, index) -> AbstractNode

Retrieve the node object at the given index within the graph.
"""
function getnode end


"""
    getedge(graph, index) -> AbstractEdge

Retrieve the edge object at the given index within the graph.
"""
function getedge(graph, index)
    throw(ErrorException("The method is not implemented for $(typeof(graph))"))
end


"""
    getedge(graph, u, v) -> AbstractEdge

Retrieve an edge object which connects the given nodes.
"""
function getedge(graph, u, v)
    throw(ErrorException("The method is not implemented for $(typeof(graph))"))
end


"""
    hasedge(graph, u, v) -> AbstractEdge

Return whether the given two nodes are connected by at least one edge.
"""
function hasedge end


"""
    nodesiter(graph)

An iterator that yields `(i, n)` where `i` is the node index, and `n` is the
node object at the index `i` within the graph.
"""
function nodesiter end


"""
    edgesiter(graph)

An iterator that yields `(i, e)` where `i` is the edge index, and `e` is the
edge object at the index `i` within the graph.
"""
function edgesiter end


"""
    nodekeys(graph) -> Vector{Int}

Return graph node keys. If the given graph is a vector graph, the keys are in
ascending order, whereas the order of indices in map graph is not guaranteed.
"""
function nodekeys end


"""
    nodeset(graph) -> Set{Int}

Return the set of node keys.
"""
function nodeset end


"""
    edgekeys(graph) -> Vector{Int}

Return graph edge keys. If the given graph is a vector graph, the keys are in
ascending order, whereas the order of indices in map graph is not guaranteed.
"""
function edgekeys end


"""
    edgeset(graph) -> Set{Int}

Return the set of edge keys.
"""
function edgeset end


"""
    neighbors(graph, n) -> Dict{Int,Int}

Return the mapping of adjacent node keys and incident edge keys connected to
the given node. If the graph is directed graph, both outneighbors and
inneighbors are mapped.
"""
function neighbors end


"""
    outneighbors(graph, n) -> Dict{Int,Int}

Return the mapping of successor node keys and out edge keys connected to
the given node.
"""
function outneighbors end


"""
    inneighbors(graph, n) -> Dict{Int,Int}

Return the mapping of predecessor node keys and in edge keys connected to
the given node.
"""
function inneighbors end


"""
    updatenode!(graph, node, n)

Rebind the node object stored at `n` of the graph by the given node object.
If the index does not exist, add new node to the position `n`.
"""
function updatenode! end


"""
    updateedge!(graph, edge, e)

Rebind the edge object stored at `e` of the graph by the given edge object.
If the index does not exist, add new edge to the position `e`.
"""
function updateedge!(graph, edge, e)
    throw(ErrorException("The method is not implemented for $(typeof(graph))"))
end


"""
    updateedge!(graph, edge, u, v)

Rebind the edge that connects nodes `u` and `v` by the given edge object.
If the nodes do not exist, throws `KeyError`.
"""
function updateedge!(graph, edge, u, v)
    throw(ErrorException("The method is not implemented for $(typeof(graph))"))
end


"""
    unlinknodes(graph, nodes)

Delete given nodes and its incident edges from the graph.
"""
function unlinknodes end


"""
    unlinkedges(graph, edges)

Delete given edges from the graph.
"""
function unlinkedges end


"""
    nodecount(graph) -> Int

Return the number of graph nodes.
"""
function nodecount end


"""
    edgecount(graph) -> Int

Return the number of graph edges.
"""
function edgecount end


"""
    neighborcount(graph, n) -> Int
    degree(graph, n) -> Int

Return the number of adjacent nodes of the node 'n'.
"""
function neighborcount end
degree = neighborcount


"""
    indegree(graph, n) -> Int

Return the number of inneighbors of the node 'n'.
"""
function indegree end


"""
    outdegree(graph, n) -> Int

Return the number of outneighbors of the node 'n'.
"""
function outdegree end


"""
    nodetype(graph)

Return the node type of the graph
"""
function nodetype end

"""
    edgetype(graph)

Return the edge type of the graph
"""
function edgetype end


# Cached properties

macro cache(ex)
    func = ex.args[1].args[1]
    dummy = gensym()
    ex.args[1].args[1] = dummy
    return quote
        $(esc(ex))
        Core.@__doc__ function $(esc(func))(graph)
            if !isdefined(graph, :cache)
                # Cache not available
                return $(esc(dummy))(graph)
            end
            symf = nameof($(esc(func)))
            if symf in keys(graph.cache)
                # Return cache
                return graph.cache[symf]
            end
            # Otherwise, set cache
            graph.cache[symf] = $(esc(dummy))(graph)
            return graph.cache[symf]
        end
    end
end


function clearcache!(graph)
    empty!(graph.cache)
    return
end
