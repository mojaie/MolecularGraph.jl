#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    PlainHyperGraph,
    plainhypergraph


struct PlainHyperGraph <: OrderedHyperGraph
    incidences::Vector{Set{Int}}
    edges::Vector{Set{Int}}
    cache::Dict{Symbol,Any}
end


"""
    plainhypergraph() -> PlainHyperGraph

Generate empty `PlainGraph`.
"""
plainhypergraph() = PlainHyperGraph([], [], Dict())

"""
    plainhypergraph(size, edges) -> PlainHyperGraph

Generate hyper graph that have given size (number of nodes) and edges (Sets of
nodes).
"""
function plainhypergraph(size::Int, edges)
    incidences = [Set{Int}() for i in 1:size]
    edges = collect(edges)
    for (i, ns) in enumerate(edges)
        for n in ns
            push!(incidences[n], i)
        end
    end
    return PlainHyperGraph(incidences, edges, Dict{Symbol,Any}())
end

"""
    plainhypergraph(graph::OrderedHyperGraph) -> PlainHyperGraph

Convert an arbitrary `OrderedHyperGraph` into a `PlainHyperGraph`.
"""
function plainhypergraph(graph::OrderedHyperGraph)
    newg = plainhypergraph()
    for incset in graph.incidences
        push!(newg.incidences, copy(incset))
    end
    for eset in graph.edges
        push!(newg.edges, copy(eset))
    end
    return newg
end

"""
    plainhypergraph(graph::HyperGraph) -> PlainHyperGraph

Convert an arbitrary `HyperGraph` into a `PlainHyperGraph`.

If the given graph has non-sequential indices (ex. subgraph view), node and
edge indices are sorted in ascending order and are re-indexed.
"""
function plainhypergraph(graph::HyperGraph)
    newg = plainhypergraph()
    nkeys = sort(collect(nodeset(graph)))
    ekeys = sort(collect(edgeset(graph)))
    nmap = Dict{Int,Int}()
    for (i, n) in enumerate(nkeys)
        nmap[n] = i
        push!(newg.incidences, Set())
    end
    for (i, e) in enumerate(ekeys)
        ns = Set{Int}()
        for oldn in getedge(graph, e)
            push!(ns, nmap[oldn])
            push!(newg.incidences[nmap[oldn]], i)
        end
        push!(newg.edges, ns)
    end
    return newg
end
