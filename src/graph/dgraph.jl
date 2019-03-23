#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    Arrow, MapDiGraph,
    setnodes, mapdigraph


struct Arrow <: DirectedEdge
    source::Int
    target::Int
end

setnodes(arrow::Arrow, s, t) = Arrow(s, t)


struct MapDiGraph{N<:AbstractNode,E<:DirectedEdge} <: DiGraph
    nodes::Dict{Int,N}
    edges::Dict{Int,E}
    outneighbormap::Dict{Int,Dict{Int,Int}}
    inneighbormap::Dict{Int,Dict{Int,Int}}

    function MapDiGraph{N,E}() where {N<:AbstractNode,E<:DirectedEdge}
        new(Dict(), Dict(), Dict(), Dict())
    end
end

"""
    mapdigraph(::Type{N}, ::Type{E}
        ) where {N<:AbstractNode,E<:DirectedEdge} -> MapDiGraph{N,E}()

Generate empty `MapDiGraph` that have nodes and edges with the given types.
"""
mapdigraph(::Type{N}, ::Type{E}
    ) where {N<:AbstractNode,E<:DirectedEdge} = MapDiGraph{N,E}()

"""
    mapdigraph(nodes, edges) -> MapDiGraph{Node,Arrow}

Generate `MapDiGraph` that have given nodes and edges represented by the list of
node indices in integer and the list of pairs of node indices, respectively.
"""
function mapdigraph(nodes, edges)
    graph = MapDiGraph{Node,Arrow}()
    for node in nodes
        graph.nodes[node] = Node()
        graph.outneighbormap[node] = Dict()
        graph.inneighbormap[node] = Dict()
    end
    for (i, edge) in enumerate(edges)
        graph.edges[i] = Arrow(edge...)
        graph.outneighbormap[edge[1]][edge[2]] = i
        graph.inneighbormap[edge[2]][edge[1]] = i
    end
    return graph
end


"""
    mapdigraph(graph::DirectedGraph; clone=false) -> MapDiGraph

Convert the given graph into a new `MapDiGraph`. The node type and edge type are
inherited from the original graph. If the given graph is `MapDiGraph`, return
a copy of the graph.


"""
function mapdigraph(graph::DirectedGraph)
    newg = mapdigraph(nodetype(graph), edgetype(graph))
    for (i, node) in nodesiter(graph)
        newg.nodes[i] = clone(node)
        newg.outneighbormap[i] = Dict()
        newg.inneighbormap[i] = Dict()
    end
    for (i, edge) in edgesiter(graph)
        newg.edges[i] = setnodes(edge, edge.source, edge.target)
        newg.outneighbormap[edge.source][edge.target] = i
        newg.inneighbormap[edge.target][edge.source] = i
    end
    return newg
end


getnode(graph::DiGraph, idx) = graph.nodes[idx]
getedge(graph::DiGraph, idx) = graph.edges[idx]
getedge(graph::DiGraph, s, t) = getedge(graph, graph.outneighbormap[s][t])
hasedge(graph::DiGraph, s, t) = haskey(graph.outneighbormap[s], t)

nodesiter(graph::MapDiGraph) = graph.nodes
edgesiter(graph::MapDiGraph) = graph.edges
nodekeys(graph::MapDiGraph) = collect(keys(graph.nodes))
edgekeys(graph::MapDiGraph) = collect(keys(graph.edges))
nodeset(graph::MapDiGraph) = Set(keys(graph.nodes))
edgeset(graph::MapDiGraph) = Set(keys(graph.edges))

outneighbors(graph::DiGraph, i) = graph.outneighbormap[i]
inneighbors(graph::DiGraph, i) = graph.inneighbormap[i]


function updatenode!(graph::MapDiGraph, node, idx)
    graph.nodes[idx] = node
    if !haskey(graph.outneighbormap, idx)
        graph.outneighbormap[idx] = Dict()
        graph.inneighbormap[idx] = Dict()
    end
    return
end

function updatenode!(graph::MapDiGraph, node)
    i = maximum(nodeset(graph)) + 1
    graph.nodes[i] = node
    graph.outneighbormap[i] = Dict()
    graph.inneighbormap[i] = Dict()
    return
end


function updateedge!(graph::MapDiGraph, edge, idx)
    nodes = nodeset(graph)
    (edge.source in nodes) || throw(KeyError(edge.source))
    (edge.target in nodes) || throw(KeyError(edge.target))
    if haskey(graph.edges, idx)
        old = getedge(graph, idx)
        delete!(graph.outneighbormap[old.source], old.target)
        delete!(graph.inneighbormap[old.target], old.source)
    end
    graph.edges[idx] = edge
    graph.outneighbormap[edge.source][edge.target] = idx
    graph.inneighbormap[edge.target][edge.source] = idx
    return
end

updateedge!(graph::MapDiGraph, edge
    ) = updateedge!(graph, edge, maximum(nodeset(graph)) + 1)

updateedge!(graph::MapDiGraph, edge, s, t
    ) = updateedge!(graph, edge, outneighbors(graph, s)[t])


function unlinknode!(graph::MapDiGraph, idx)
    (idx in nodeset(graph)) || throw(KeyError(idx))
    for (s, succ) in outneighbors(graph, idx)
        delete!(graph.edges, succ)
        delete!(graph.inneighbormap[s], idx)
    end
    for (p, pred) in inneighbors(graph, idx)
        delete!(graph.edges, pred)
        delete!(graph.outneighbormap[p], idx)
    end
    delete!(graph.nodes, idx)
    delete!(graph.outneighbormap, idx)
    delete!(graph.inneighbormap, idx)
    return
end


function unlinkedge!(graph::MapDiGraph, source, target)
    ns = nodeset(graph)
    (source in ns) || throw(KeyError(source))
    (target in ns) || throw(KeyError(target))
    delete!(graph.edges, graph.outneighbormap[source][target])
    delete!(graph.outneighbormap[source], target)
    delete!(graph.inneighbormap[target], source)
    return
end

function unlinkedge!(graph::MapDiGraph, idx)
    (idx in nodeset(graph)) || throw(KeyError(idx))
    a = getedge(graph, idx)
    delete!(graph.edges, idx)
    delete!(graph.outneighbormap[a.source], a.target)
    delete!(graph.inneighbormap[a.target], a.source)
    return
end


nodetype(graph::MapDiGraph) = valtype(graph.nodes)
edgetype(graph::MapDiGraph) = valtype(graph.edges)
