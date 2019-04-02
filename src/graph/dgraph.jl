#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    Arrow, DiGraph,
    setnodes, digraph


struct Arrow <: DirectedEdge
    source::Int
    target::Int
end

setnodes(arrow::Arrow, s, t) = Arrow(s, t)


struct DiGraph{N<:AbstractNode,E<:DirectedEdge} <: DirectedGraph
    nodes::Vector{N}
    edges::Vector{E}
    outneighbormap::Vector{Dict{Int,Int}}
    inneighbormap::Vector{Dict{Int,Int}}

    function DiGraph{N,E}() where {N<:AbstractNode,E<:DirectedEdge}
        new([], [], [], [])
    end
end

"""
    digraph(::Type{N}, ::Type{E}
        ) where {N<:AbstractNode,E<:DirectedEdge} -> DiGraph{N,E}()

Generate empty `DiGraph` that have nodes and edges with the given types.
"""
digraph(::Type{N}, ::Type{E}
    ) where {N<:AbstractNode,E<:DirectedEdge} = DiGraph{N,E}()

"""
    digraph(nodes, edges) -> DiGraph{Node,Arrow}

Generate `DiGraph` that have given nodes and edges represented by the list of
node indices in integer and the list of pairs of node indices, respectively.
"""
function digraph(size::Int, edges)
    graph = DiGraph{Node,Arrow}()
    append!(graph.nodes, [Node() for i in 1:size])
    append!(graph.inneighbormap, [Dict() for i in 1:size])
    append!(graph.outneighbormap, [Dict() for i in 1:size])
    for (i, (s, t)) in enumerate(edges)
        push!(graph.edges, Arrow(s, t))
        graph.outneighbormap[s][t] = i
        graph.inneighbormap[t][s] = i
    end
    return graph
end


"""
    digraph(graph::DirectedGraph) -> DiGraph

Convert the given graph into a new `DiGraph`. The node type and edge type are
inherited from the original graph. If the given graph has non-sequential
indices (ex. subgraph view), node and edge indices are sorted in ascending
order and are re-indexed.

If you really need mutable nodes and edges, and want the graph with fully deepcopied elements, implement `clone` and `setnodes` as deepcopy methods for
nodes and edges, respectively.
"""
function digraph(graph::DirectedGraph)
    newg = digraph(nodetype(graph), edgetype(graph))
    nkeys = sort(nodekeys(graph))
    ekeys = sort(edgekeys(graph))
    nmap = Dict{Int,Int}()
    for (i, n) in enumerate(nkeys)
        nmap[n] = i
        node = getnode(graph, n)
        push!(newg.nodes, clone(node))
        push!(newg.outneighbormap, Dict())
        push!(newg.inneighbormap, Dict())
    end
    for (i, e) in enumerate(ekeys)
        edge = getedge(graph, e)
        s = nmap[edge.source]
        t = nmap[edge.target]
        push!(newg.edges, setnodes(edge, s, t))
        newg.outneighbormap[s][t] = i
        newg.inneighbormap[t][s] = i
    end
    return newg
end


getnode(graph::DiGraph, idx) = graph.nodes[idx]
getedge(graph::DiGraph, idx) = graph.edges[idx]
getedge(graph::DiGraph, s, t) = getedge(graph, graph.outneighbormap[s][t])
hasedge(graph::DiGraph, s, t) = haskey(graph.outneighbormap[s], t)
outneighbors(graph::DiGraph, i) = graph.outneighbormap[i]
inneighbors(graph::DiGraph, i) = graph.inneighbormap[i]

nodekeys(graph::DiGraph) = collect(1:nodecount(graph))
edgekeys(graph::DiGraph) = collect(1:edgecount(graph))
nodevalues(graph::DiGraph) = graph.nodes
edgevalues(graph::DiGraph) = graph.edges
nodesiter(graph::DiGraph) = enumerate(graph.nodes)
edgesiter(graph::DiGraph) = enumerate(graph.edges)
nodeset(graph::DiGraph) = Set(1:nodecount(graph))
edgeset(graph::DiGraph) = Set(1:edgecount(graph))

nodecount(graph::DiGraph) = length(graph.nodes)
edgecount(graph::DiGraph) = length(graph.edges)


function updatenode!(graph::DiGraph, node, idx)
    idx > nodecount(graph) + 1 && throw(DomainError(idx))
    if idx == nodecount(graph) + 1
        updatenode!(graph, node)
        return
    end
    graph.nodes[idx] = node
    return
end

function updatenode!(graph::DiGraph, node)
    push!(graph.nodes, node)
    push!(graph.outneighbormap, Dict())
    push!(graph.inneighbormap, Dict())
    return
end


function updateedge!(graph::DiGraph, edge, idx)
    idx > edgecount(graph) + 1 && throw(DomainError(idx))
    if idx == edgecount(graph) + 1
        updateedge!(graph, edge)
        return
    end
    ncnt = nodecount(graph)
    edge.source > ncnt && throw(DomainError(edge.source))
    edge.target > ncnt && throw(DomainError(edge.target))
    old = getedge(graph, idx)
    delete!(graph.outneighbormap[old.source], old.target)
    delete!(graph.inneighbormap[old.target], old.source)
    graph.edges[idx] = edge
    graph.outneighbormap[edge.source][edge.target] = idx
    graph.inneighbormap[edge.target][edge.source] = idx
    return
end

function updateedge!(graph::DiGraph, edge)
    ncnt = nodecount(graph)
    edge.source > ncnt && throw(DomainError(edge.source))
    edge.target > ncnt && throw(DomainError(edge.target))
    i = edgecount(graph) + 1
    push!(graph.edges, edge)
    graph.outneighbormap[edge.source][edge.target] = i
    graph.inneighbormap[edge.target][edge.source] = i
    return
end

updateedge!(graph::DiGraph, edge, s, t
    ) = updateedge!(graph, edge, outneighbors(graph, s)[t])


function unlinknodes(graph::DiGraph, nodes)
    # TODO: costful
    subg = nodesubgraph(graph, setdiff(nodeset(graph), nodes))
    return digraph(subg)
end


function unlinkedges(graph::DiGraph, edges)
    # TODO: costful
    subg = DiSubgraphView(graph, nodeset(graph), setdiff(edgeset(graph), edges))
    return digraph(subg)
end


nodetype(graph::DiGraph) = eltype(graph.nodes)
edgetype(graph::DiGraph) = eltype(graph.edges)
