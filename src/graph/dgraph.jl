#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    Arrow,
    MapDGraph,
    connect,
    getnode,
    getedge,
    nodesiter,
    edgesiter,
    nodekeys,
    edgekeys,
    predecessors,
    successors,
    nodecount,
    edgecount,
    succkeys,
    predkeys,
    succnodes,
    prednodes,
    inedgekeys,
    outedgekeys,
    inedges,
    outedges,
    indegree,
    outdegree,
    updatenode!,
    updateedge!,
    unlinknode!,
    unlinkedge!,
    similarmap



struct Arrow <: AbstractDirectedEdge
    source::Int
    target::Int
    attr::Dict
end

Arrow(s, t) = Arrow(s, t, Dict())
connect(a::Arrow, s, t) = Arrow(s, t, a.attr)


struct MapDGraph{N<:AbstractNode,E<:AbstractDirectedEdge} <: DGraph
    nodes::Dict{Int,N}
    edges::Dict{Int,E}
    successors::Dict{Int,Dict{Int,Int}}
    predecessors::Dict{Int,Dict{Int,Int}}

    function MapDGraph{N,E}() where {N<:AbstractNode,E<:AbstractDirectedEdge}
        new(Dict(), Dict(), Dict(), Dict())
    end
end

function MapDGraph(nodes::AbstractArray{Int},
                    edges::AbstractArray{Tuple{Int,Int}})
    graph = MapDGraph{Node,Arrow}()
    for node in nodes
        updatenode!(graph, Node(), node)
    end
    for (i, edge) in enumerate(edges)
        updateedge!(graph, Arrow(edge...), i)
    end
    graph
end

getnode(graph::MapDGraph, idx) = graph.nodes[idx]

getedge(graph::MapDGraph, idx) = graph.edges[idx]
getedge(graph::MapDGraph, s, t) = getedge(graph, graph.successors[s][t])

nodesiter(graph::MapDGraph) = graph.nodes

nodekeys(graph::MapDGraph) = Set(keys(graph.nodes))

edgesiter(graph::MapDGraph) = graph.edges

edgekeys(graph::MapDGraph) = Set(keys(graph.edges))

successors(graph::MapDGraph, idx) = graph.successors[idx]

predecessors(graph::MapDGraph, idx) = graph.predecessors[idx]

nodecount(graph::MapDGraph) = length(graph.nodes)

edgecount(graph::MapDGraph) = length(graph.edges)

succkeys(graph::AbstractDGraph, idx) = collect(keys(successors(graph, idx)))
predkeys(graph::AbstractDGraph, idx) = collect(keys(predecessors(graph, idx)))

succnodes(
    graph::AbstractDGraph, idx) = getnode.((graph,), succkeys(graph, idx))
prednodes(
    graph::AbstractDGraph, idx) = getnode.((graph,), predkeys(graph, idx))

outedgekeys(
    graph::AbstractDGraph, idx) = collect(values(successors(graph, idx)))
inedgekeys(
    graph::AbstractDGraph, idx) = collect(values(predecessors(graph, idx)))

outedges(
    graph::AbstractDGraph, idx) = getedge.((graph,), outedgekeys(graph, idx))
inedges(
    graph::AbstractDGraph, idx) = getedge.((graph,), inedgekeys(graph, idx))

outdegree(graph::AbstractDGraph, idx) = length(successors(graph, idx))
indegree(graph::AbstractDGraph, idx) = length(predecessors(graph, idx))


function updatenode!(graph::MapDGraph, node, idx)
    """Add or update a node"""
    graph.nodes[idx] = node
    if !(idx in keys(graph.successors))
        graph.successors[idx] = Dict{Int,Int}()
        graph.predecessors[idx] = Dict{Int,Int}()
    end
    return
end


function updateedge!(graph::MapDGraph, edge, idx)
    """Add or update an edge"""
    if !(edge.source in keys(graph.nodes))
        throw(OperationError("Missing node: $(edge.source)"))
    elseif !(edge.target in keys(graph.nodes))
        throw(OperationError("Missing node: $(edge.target)"))
    end
    graph.edges[idx] = edge
    graph.successors[edge.source][edge.target] = idx
    graph.predecessors[edge.target][edge.source] = idx
    return
end

updateedge!(
    G::MapDGraph, edge, s, t) = updateedge!(G, edge, graph.successors[s][t])


function unlinknode!(graph::MapDGraph, idx)
    """Remove a node and its connecting edges"""
    if !(idx in keys(graph.nodes))
        throw(OperationError("Missing node: $(idx)"))
    end
    for (s, succ) in graph.successors[idx]
        delete!(graph.edges, succ)
        delete!(graph.predecessors[s], idx)
    end
    for (p, pred) in graph.predecessors[idx]
        delete!(graph.edges, pred)
        delete!(graph.successors[p], idx)
    end
    delete!(graph.nodes, idx)
    delete!(graph.successors, idx)
    delete!(graph.predecessors, idx)
    return
end


function unlinkedge!(graph::MapDGraph, source, target)
    """Remove an edge"""
    if !(source in keys(graph.nodes))
        throw(OperationError("Missing node: $(source)"))
    elseif !(target in keys(graph.nodes))
        throw(OperationError("Missing node: $(target)"))
    end
    delete!(graph.edges, graph.successors[source][target])
    delete!(graph.successors[source], target)
    delete!(graph.predecessors[target], source)
    return
end

function unlinkedge!(graph::MapDGraph, idx)
    """Remove an edge"""
    if !(idx in keys(graph.edges))
        throw(OperationError("Missing edge: $(idx)"))
    end
    a = getedge(graph, idx)
    delete!(graph.edges, idx)
    delete!(graph.successors[a.source], a.target)
    delete!(graph.predecessors[a.target], a.source)
    return
end


function similarmap(graph::MapDGraph)
    N = valtype(graph.nodes)
    E = valtype(graph.edges)
    MapDGraph{N,E}()
end
