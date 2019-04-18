#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    PlainDiGraph, plaindigraph


struct PlainDiGraph <: OrderedDiGraph
    outneighbormap::Vector{Dict{Int,Int}}
    inneighbormap::Vector{Dict{Int,Int}}
    edges::Vector{Tuple{Int,Int}}
    cache::Dict{Symbol,Any}
end

"""
    plaindigraph() -> PlainDiGraph

Generate empty `PlainDiGraph`.
"""
plaindigraph() = PlainDiGraph([], [], [], Dict())

"""
    plaindigraph(nodes, edges) -> PlainDiGraph

Generate `DiGraph` that have given nodes and edges represented by the list of
node indices in integer and the list of pairs of node indices, respectively.
"""
function plaindigraph(size::Int, edges)
    outnbrmap = [Dict{Int,Int}() for i in 1:size]
    innbrmap = [Dict{Int,Int}() for i in 1:size]
    edges = collect(edges)
    for (i, (s, t)) in enumerate(edges)
        outnbrmap[s][i] = t
        innbrmap[t][i] = s
    end
    return PlainDiGraph(outnbrmap, innbrmap, edges, Dict{Symbol,Any}())
end

"""
    plaindigraph(graph::OrderedDiGraph) -> PlainDiGraph

Convert an arbitrary `OrderedDiGraph` into a `PlainDiGraph`.
"""
function plaindigraph(graph::OrderedDiGraph)
    newg = plaindigraph()
    for nbr in graph.outneighbormap
        push!(newg.outneighbormap, copy(nbr))
    end
    for nbr in graph.inneighbormap
        push!(newg.inneighbormap, copy(nbr))
    end
    append!(newg.edges, graph.edges)
    return newg
end

"""
    plaindigraph(graph::DirectedGraph) -> PlainDiGraph

Convert an arbitrary `DirectedGraph` into a `PlainDiGraph`.

If the given graph has non-sequential indices (ex. subgraph view), node and
edge indices are sorted in ascending order and are re-indexed.
"""
function plaindigraph(graph::DirectedGraph)
    newg = plaindigraph()
    nkeys = sort(collect(nodeset(graph)))
    ekeys = sort(collect(edgeset(graph)))
    nmap = Dict{Int,Int}()
    for (i, n) in enumerate(nkeys)
        nmap[n] = i
        push!(newg.outneighbormap, Dict())
        push!(newg.inneighbormap, Dict())
    end
    for (i, e) in enumerate(ekeys)
        (olds, oldt) = getedge(graph, e)
        s = nmap[olds]
        t = nmap[oldt]
        push!(newg.edges, (s, t))
        newg.outneighbormap[s][i] = t
        newg.inneighbormap[t][i] = s
    end
    return newg
end


function unlinknodes(graph::PlainDiGraph, nodes)
    # TODO: costful
    subg = nodesubgraph(graph, setdiff(nodeset(graph), nodes))
    return plaindigraph(subg)
end


function unlinkedges(graph::PlainDiGraph, edges)
    # TODO: costful
    subg = SubgraphView(graph, nodeset(graph), setdiff(edgeset(graph), edges))
    return plaindigraph(subg)
end
