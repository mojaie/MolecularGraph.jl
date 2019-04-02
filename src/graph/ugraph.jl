#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    Edge, VectorGraph, MapGraph,
    setnodes, vectorgraph, mapgraph


struct Edge <: UndirectedEdge
    u::Int
    v::Int
end

setnodes(edge::Edge, u, v) = Edge(u, v)


struct VectorGraph{N<:AbstractNode,E<:UndirectedEdge} <: UndirectedGraph
    nodes::Vector{N}
    edges::Vector{E}
    neighbormap::Vector{Dict{Int,Int}}
    cache::Dict{Symbol,Any}

    function VectorGraph{N,E}() where {N<:AbstractNode,E<:UndirectedEdge}
        new([], [], [], Dict())
    end
end


struct MapGraph{N<:AbstractNode,E<:UndirectedEdge} <: UndirectedGraph
    # TODO: Deprecated
    nodes::Dict{Int,N}
    edges::Dict{Int,E}
    neighbormap::Dict{Int,Dict{Int,Int}}

    function MapGraph{N,E}() where {N<:AbstractNode,E<:UndirectedEdge}
        new(Dict(), Dict(), Dict())
    end
end


UGraph = Union{VectorGraph,MapGraph}


"""
    vectorgraph(::Type{N}, ::Type{E}
        ) where {N<:AbstractNode,E<:UndirectedEdge} -> VectorGraph{N,E}()

Generate empty vector graph that have nodes and edges with the given types.
"""
vectorgraph(::Type{N}, ::Type{E}
    ) where {N<:AbstractNode,E<:UndirectedEdge} = VectorGraph{N,E}()

"""
    vectorgraph(nodes, edges) -> VectorGraph{Node,Edge}

Generate vector graph that have given nodes and edges represented by the list of
node indices in integer and the list of pairs of node indices, respectively.
"""
function vectorgraph(size::Int, edges)
    graph = VectorGraph{Node,Edge}()
    append!(graph.nodes, [Node() for i in 1:size])
    append!(graph.neighbormap, [Dict() for i in 1:size])
    for (i, (u, v)) in enumerate(edges)
        push!(graph.edges, Edge(u, v))
        graph.neighbormap[u][v] = i
        graph.neighbormap[v][u] = i
    end
    return graph
end

"""
    vectorgraph(graph::UndirectedGraph) -> VectorGraph

Convert the given graph into a new `VectorGraph`. The node type and edge type
are inherited from the original graph. If the given graph has non-sequential
indices (ex. subgraph view), node and edge indices are sorted in ascending
order and are re-indexed.

If you really need mutable nodes and edges, and want the graph with fully
deepcopied elements, implement `clone` and `setnodes` as deepcopy methods for
nodes and edges, respectively.
"""
function vectorgraph(graph::UndirectedGraph)
    newg = vectorgraph(nodetype(graph), edgetype(graph))
    nkeys = sort(nodekeys(graph))
    ekeys = sort(edgekeys(graph))
    nmap = Dict{Int,Int}()
    for (i, n) in enumerate(nkeys)
        nmap[n] = i
        node = getnode(graph, n)
        push!(newg.nodes, clone(node))
        push!(newg.neighbormap, Dict())
    end
    for (i, e) in enumerate(ekeys)
        edge = getedge(graph, e)
        u = nmap[edge.u]
        v = nmap[edge.v]
        push!(newg.edges, setnodes(edge, u, v))
        newg.neighbormap[u][v] = i
        newg.neighbormap[v][u] = i
    end
    return newg
end


Base.getindex(graph::VectorGraph, sym::Symbol) = eval(Expr(:call, sym, graph))
Base.getindex(
    graph::VectorGraph, k1::Symbol, k2::Symbol, K::Symbol...
) = hcat(eval(Expr(:call, sym, graph)) for k in [k1, k2, K...])

getnode(graph::UGraph, idx) = graph.nodes[idx]
getedge(graph::UGraph, idx) = graph.edges[idx]
getedge(graph::UGraph, u, v) = getedge(graph, graph.neighbormap[u][v])
hasedge(graph::UGraph, u, v) = haskey(graph.neighbormap[u], v)
neighbors(graph::UGraph, idx) = graph.neighbormap[idx]

nodekeys(graph::VectorGraph) = collect(1:nodecount(graph))
edgekeys(graph::VectorGraph) = collect(1:edgecount(graph))
nodevalues(graph::VectorGraph) = graph.nodes
edgevalues(graph::VectorGraph) = graph.edges
nodesiter(graph::VectorGraph) = enumerate(graph.nodes)
edgesiter(graph::VectorGraph) = enumerate(graph.edges)
nodeset(graph::VectorGraph) = Set(1:nodecount(graph))
edgeset(graph::VectorGraph) = Set(1:edgecount(graph))

nodecount(graph::UGraph) = length(graph.nodes)
edgecount(graph::UGraph) = length(graph.edges)


function updatenode!(graph::VectorGraph, node, idx)
    idx > nodecount(graph) + 1 && throw(DomainError(idx))
    if idx == nodecount(graph) + 1
        updatenode!(graph, node)
        return
    end
    graph.nodes[idx] = node
    return
end

function updatenode!(graph::VectorGraph, node)
    push!(graph.nodes, node)
    push!(graph.neighbormap, Dict())
    return
end


function updateedge!(graph::VectorGraph, edge, idx)
    idx > edgecount(graph) + 1 && throw(DomainError(idx))
    if idx == edgecount(graph) + 1
        updateedge!(graph, edge)
        return
    end
    ncnt = nodecount(graph)
    edge.u > ncnt && throw(DomainError(edge.u))
    edge.v > ncnt && throw(DomainError(edge.v))
    old = getedge(graph, idx)
    delete!(graph.neighbormap[old.u], old.v)
    delete!(graph.neighbormap[old.v], old.u)
    graph.edges[idx] = edge
    graph.neighbormap[edge.u][edge.v] = idx
    graph.neighbormap[edge.v][edge.u] = idx
    return
end

function updateedge!(graph::VectorGraph, edge)
    ncnt = nodecount(graph)
    edge.u > ncnt && throw(DomainError(edge.u))
    edge.v > ncnt && throw(DomainError(edge.v))
    i = edgecount(graph) + 1
    push!(graph.edges, edge)
    graph.neighbormap[edge.u][edge.v] = i
    graph.neighbormap[edge.v][edge.u] = i
    return
end

updateedge!(graph::UGraph, edge, u, v
    ) = updateedge!(graph, edge, neighbors(graph, u)[v])


function unlinknodes(graph::VectorGraph, nodes)
    # TODO: costful
    subg = nodesubgraph(graph, setdiff(nodeset(graph), nodes))
    return vectorgraph(subg)
end


function unlinkedges(graph::VectorGraph, edges)
    # TODO: costful
    subg = SubgraphView(graph, nodeset(graph), setdiff(edgeset(graph), edges))
    return vectorgraph(subg)
end


nodetype(graph::VectorGraph) = eltype(graph.nodes)
edgetype(graph::VectorGraph) = eltype(graph.edges)




mapgraph(::Type{N}, ::Type{E}
    ) where {N<:AbstractNode,E<:UndirectedEdge} = MapGraph{N,E}()

function mapgraph(nodes, edges)
    graph = MapGraph{Node,Edge}()
    for node in nodes
        graph.nodes[node] = Node()
        graph.neighbormap[node] = Dict()
    end
    for (i, edge) in enumerate(edges)
        (u, v) = edge
        graph.edges[i] = Edge(u, v)
        graph.neighbormap[u][v] = i
        graph.neighbormap[v][u] = i
    end
    return graph
end

function mapgraph(graph::UndirectedGraph)
    newg = mapgraph(nodetype(graph), edgetype(graph))
    for (i, node) in nodesiter(graph)
        newg.nodes[i] = clone(node)
        newg.neighbormap[i] = Dict()
    end
    for (i, edge) in edgesiter(graph)
        newg.edges[i] = setnodes(edge, edge.u, edge.v)
        newg.neighbormap[edge.u][edge.v] = i
        newg.neighbormap[edge.v][edge.u] = i
    end
    return newg
end


nodekeys(graph::MapGraph) = collect(keys(graph.nodes))
edgekeys(graph::MapGraph) = collect(keys(graph.edges))
# TODO: `enumerate` yields `Tuple` whereas `Dict` yields `Pair`
nodesiter(graph::MapGraph) = graph.nodes
edgesiter(graph::MapGraph) = graph.edges
nodeset(graph::MapGraph) = Set(keys(graph.nodes))
edgeset(graph::MapGraph) = Set(keys(graph.edges))

function updatenode!(graph::MapGraph, node, idx)
    graph.nodes[idx] = node
    if !haskey(graph.neighbormap, idx)
        graph.neighbormap[idx] = Dict()
    end
    return
end

function updatenode!(graph::MapGraph, node)
    i = maximum(nodeset(graph)) + 1
    graph.nodes[i] = node
    graph.neighbormap[i] = Dict()
    return
end

function updateedge!(graph::MapGraph, edge, idx)
    nodes = nodeset(graph)
    (edge.u in nodes) || throw(KeyError(edge.u))
    (edge.v in nodes) || throw(KeyError(edge.v))
    if haskey(graph.edges, idx)
        old = getedge(graph, idx)
        delete!(graph.neighbormap[old.u], old.v)
        delete!(graph.neighbormap[old.v], old.u)
    end
    graph.edges[idx] = edge
    graph.neighbormap[edge.u][edge.v] = idx
    graph.neighbormap[edge.v][edge.u] = idx
    return
end

updateedge!(graph::MapGraph, edge
    ) = updateedge!(graph, edge, maximum(edgeset(graph)) + 1)

"""
function unlinknode!(graph::MapGraph, idx)
    (idx in nodekeys(graph)) || throw(KeyError(idx))
    for (n, nbr) in neighbors(graph, idx)
        delete!(graph.edges, nbr)
        delete!(graph.neighbormap[n], idx)
    end
    delete!(graph.nodes, idx)
    delete!(graph.neighbormap, idx)
    return
end

function unlinkedge!(graph::MapGraph, u, v)
    nodes = nodekeys(graph)
    (u in nodes) || throw(KeyError(u))
    (v in nodes) || throw(KeyError(v))
    delete!(graph.edges, graph.neighbormap[u][v])
    delete!(graph.neighbormap[u], v)
    delete!(graph.neighbormap[v], u)
    return
end

function unlinkedge!(graph::MapGraph, idx)
    (idx in keys(graph.edges)) || throw(KeyError(idx))
    e = getedge(graph, idx)
    delete!(graph.edges, idx)
    delete!(graph.neighbormap[e.u], e.v)
    delete!(graph.neighbormap[e.v], e.u)
    return
end
"""
nodetype(graph::MapGraph) = valtype(graph.nodes)
edgetype(graph::MapGraph) = valtype(graph.edges)
