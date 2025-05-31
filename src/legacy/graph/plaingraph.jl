#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    PlainGraph, ImmutablePlainGraph,
    plaingraph, immutableplaingraph


struct PlainGraph <: OrderedGraph
    neighbormap::Vector{Dict{Int,Int}}
    edges::Vector{Tuple{Int,Int}}
    cache::Dict{Symbol,Any}
end


struct ImmutablePlainGraph <: OrderedGraph
    neighbormap::Vector{Dict{Int,Int}}
    edges::Vector{Tuple{Int,Int}}
    cache::Dict{Symbol,Any}
end


"""
    plaingraph() -> PlainGraph

Generate empty `PlainGraph`.
"""
plaingraph() = PlainGraph([], [], Dict())

"""
    plaingraph(nodes, edges) -> PlainGraph{PlainNode,Edge}

Generate vector graph that have given nodes and edges represented by the list of
node indices in integer and the list of pairs of node indices, respectively.
"""
function plaingraph(size::Int, edges)
    nbrmap = [Dict{Int,Int}() for i in 1:size]
    edges = collect(edges)
    for (i, (u, v)) in enumerate(edges)
        nbrmap[u][i] = v
        nbrmap[v][i] = u
    end
    return PlainGraph(nbrmap, edges, Dict{Symbol,Any}())
end

"""
    plaingraph(graph::OrderedGraph) -> PlainGraph

Convert an arbitrary `OrderedGraph` into a `PlainGraph`.
"""
function plaingraph(graph::OrderedGraph)
    newg = plaingraph()
    for nbr in graph.neighbormap
        push!(newg.neighbormap, copy(nbr))
    end
    append!(newg.edges, graph.edges)
    return newg
end

"""
    plaingraph(graph::UndirectedGraph) -> PlainGraph

Convert an arbitrary `UndirectedGraph` into a `PlainGraph`.

If the given graph has non-sequential indices (ex. subgraph view), node and
edge indices are sorted in ascending order and are re-indexed.
"""
function plaingraph(graph::UndirectedGraph)
    newg = plaingraph()
    nkeys = sort(collect(nodeset(graph)))
    ekeys = sort(collect(edgeset(graph)))
    nmap = Dict{Int,Int}()
    for (i, n) in enumerate(nkeys)
        nmap[n] = i
        push!(newg.neighbormap, Dict())
    end
    for (i, e) in enumerate(ekeys)
        (oldu, oldv) = getedge(graph, e)
        u = nmap[oldu]
        v = nmap[oldv]
        push!(newg.edges, (u, v))
        newg.neighbormap[u][i] = v
        newg.neighbormap[v][i] = u
    end
    return newg
end


"""
    immutableplaingraph() -> ImmutablePlainGraph
    immutableplaingraph(size::Int, edges) -> ImmutablePlainGraph
    immutableplaingraph(graph::UndirectedGraph) -> ImmutablePlainGraph

Convert an arbitrary `UndirectedGraph` into a `ImmutablePlainGraph`.

`ImmutablePlainGraph` is same as `PlainGraph` except for availability of
`addnode!` and `addedge!`, and no performance difference exists between these
graphs. `ImmutablePlainGraph` is mainly used for unit testing to check adverse
effect of some graph algorithm implementations.
"""
immutableplaingraph() = ImmutablePlainGraph([], [], Dict())
function immutableplaingraph(size::Int, edges)
    g = plaingraph(size, edges)
    return ImmutablePlainGraph(g.neighbormap, g.edges, g.cache)
end
function immutableplaingraph(graph::UndirectedGraph)
    g = plaingraph(graph)
    return ImmutablePlainGraph(g.neighbormap, g.edges, g.cache)
end


function addnode!(graph::ImmutablePlainGraph)
    throw(ErrorException("Immutable"))
end

function addedge!(graph::ImmutablePlainGraph, u::Int, v::Int)
    throw(ErrorException("Immutable"))
end


function unlinknodes(graph::PlainGraph, nodes)
    # TODO: costful
    subg = nodesubgraph(graph, setdiff(nodeset(graph), nodes))
    return plaingraph(subg)
end


function unlinkedges(graph::PlainGraph, edges)
    # TODO: costful
    subg = SubgraphView(graph, nodeset(graph), setdiff(edgeset(graph), edges))
    return plaingraph(subg)
end
