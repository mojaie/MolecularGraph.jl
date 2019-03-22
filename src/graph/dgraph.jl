#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    Arrow, similaredge, MapDGraph,
    similargraph



struct Arrow <: DirectedEdge
    source::Int
    target::Int
end

similaredge(e::Arrow, s, t) = Arrow(s, t)


struct MapDGraph{N<:AbstractNode,E<:DirectedEdge} <: DirectedGraph
    nodes::Dict{Int,N}
    edges::Dict{Int,E}
    successors::Dict{Int,Dict{Int,Int}}
    predecessors::Dict{Int,Dict{Int,Int}}

    function MapDGraph{N,E}() where {N<:AbstractNode,E<:DirectedEdge}
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
    return graph
end

getnode(graph::DirectedGraph, idx) = graph.nodes[idx]

getedge(graph::DirectedGraph, idx) = graph.edges[idx]
getedge(graph::DirectedGraph, s, t) = getedge(graph, graph.successors[s][t])

nodesiter(graph::MapDGraph) = graph.nodes
edgesiter(graph::MapDGraph) = graph.edges

nodekeys(g::MapDGraph) = collect(keys(g.nodes))
edgekeys(g::MapDGraph) = collect(keys(g.edges))

nodeset(g::MapDGraph) = Set(keys(g.nodes))
edgeset(g::MapDGraph) = Set(keys(g.edges))

successors(g::DirectedGraph, i) = g.successors[i]
predecessors(g::DirectedGraph, i) = g.predecessors[i]


function updatenode!(graph::MapDGraph, node, idx)
    graph.nodes[idx] = node
    if !(idx in keys(graph.successors))
        graph.successors[idx] = Dict{Int,Int}()
        graph.predecessors[idx] = Dict{Int,Int}()
    end
    return
end


function updateedge!(graph::MapDGraph, edge, idx)
    ns = nodekeys(graph)
    (edge.source in ns) || throw(KeyError(edge.source))
    (edge.target in ns) || throw(KeyError(edge.target))
    graph.edges[idx] = edge
    graph.successors[edge.source][edge.target] = idx
    graph.predecessors[edge.target][edge.source] = idx
    return
end

updateedge!(
    g::MapDGraph, edge, s, t) = updateedge!(g, edge, successors(graph, s)[t])


function unlinknode!(graph::MapDGraph, idx)
    (idx in keys(graph.nodes)) || throw(KeyError(idx))
    for (s, succ) in successors(graph, idx)
        delete!(graph.edges, succ)
        delete!(graph.predecessors[s], idx)
    end
    for (p, pred) in predecessors(graph, idx)
        delete!(graph.edges, pred)
        delete!(graph.successors[p], idx)
    end
    delete!(graph.nodes, idx)
    delete!(graph.successors, idx)
    delete!(graph.predecessors, idx)
    return
end


function unlinkedge!(graph::MapDGraph, source, target)
    ns = nodekeys(graph)
    (source in ns) || throw(KeyError(source))
    (target in ns) || throw(KeyError(target))
    delete!(graph.edges, graph.successors[source][target])
    delete!(graph.successors[source], target)
    delete!(graph.predecessors[target], source)
    return
end

function unlinkedge!(graph::MapDGraph, idx)
    (idx in keys(graph.edges)) || throw(KeyError(idx))
    a = getedge(graph, idx)
    delete!(graph.edges, idx)
    delete!(graph.successors[a.source], a.target)
    delete!(graph.predecessors[a.target], a.source)
    return
end


nodetype(graph::MapDGraph) = valtype(graph.nodes)
edgetype(graph::MapDGraph) = valtype(graph.edges)

function similargraph(graph::MapDGraph)
    N = nodetype(graph)
    E = edgetype(graph)
    MapDGraph{N,E}()
end
