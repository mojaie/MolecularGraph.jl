#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    getsourcenode, getsourceedge, getunionnode, getunionedge,
    disjointunion, disjointunion!


struct DisjointUnionNode <: AbstractNode
    source::Int
    sourcekey::Int
end


struct DisjointUnionEdge <: UndirectedEdge
    source::Int
    sourcekey::Int
end


struct DisjointUnionGraph <: OrderedGraph
    neighbormap::Vector{Dict{Int,Int}}
    edges::Vector{Tuple{Int,Int}}
    nodeattrs::Vector{DisjointUnionNode}
    edgeattrs::Vector{DisjointUnionEdge}
    cache::Dict{Symbol,Any}
end


getsourcenode(graph, i) = graph.nodeattrs[i]
getsourceedge(graph, i) = graph.edgeattrs[i]

function getunionnode(graph::DisjointUnionGraph, source, key)
    for (i, node) in enumerate(nodeattrs(graph))
        if node.source == source && node.sourcekey == key
            return i
        end
    end
end

function getunionedge(graph::DisjointUnionGraph, source, key)
    for (i, edge) in enumerate(edgeattrs(graph))
        if edge.source == source && edge.sourcekey == key
            return i
        end
    end
end


"""
    disjointunion(g1::UndirectedGraph, g2::UndirectedGraph,
        G::UndirectedGraph...) -> DisjointUnionGraph

Generate disjoint union graph of given graphs. The new graph with type
`DisjointUnionGraph` retains mapping to the original graphs as nodes and edges
attributes.
"""
function disjointunion(
            g1::UndirectedGraph, g2::UndirectedGraph, G::UndirectedGraph...)
    U = DisjointUnionGraph([], [], [], [], Dict())
    for (src, g) in enumerate([g1, g2, G...])
        nsort = sort(collect(nodeset(g)))
        esort = sort(collect(edgeset(g)))
        nmap = Dict(n => i + nodecount(U) for (i, n) in enumerate(nsort))
        emap = Dict(e => i + edgecount(U) for (i, e) in enumerate(esort))
        nbrs = [neighbors(g, n) for n in nsort]
        edges = [getedge(g, e) for e in esort]
        append!(U.neighbormap, relabel(nbrs, nmap, emap))
        append!(U.edges, relabel(edges, nmap))
        append!(U.nodeattrs, [DisjointUnionNode(src, k) for k in nsort])
        append!(U.edgeattrs, [DisjointUnionEdge(src, k) for k in esort])
    end
    return U
end


"""
    disjointunion!(g1::T, g2::T, G::T...) where {T<:OrderedGraph} -> T

Generate disjoint union graph of given graphs. `g1` will be overwritten by the
union graph. Unlike non-destructive `disjointunion`, `g1` does not retain any
information about other given graphs but a bit faster.
"""
function disjointunion!(g1::T, g2::T, G::T...) where {T<:OrderedGraph}
    for g in [g2, G...]
        noff = nodecount(g1)
        nmap = Dict(i => i + noff for i in 1:nodecount(g))
        eoff = edgecount(g1)
        emap = Dict(i => i + eoff for i in 1:edgecount(g))
        append!(g1.neighbormap, relabel(g.neighbormap, nmap, emap))
        append!(g1.edges, relabel(g.edges, nmap))
        isdefined(g1, :nodeattrs) && append!(g1.nodeattrs, g.nodeattrs)
        isdefined(g1, :edgeattrs) && append!(g1.edgeattrs, g.edgeattrs)
    end
end


function relabel(edges::Vector{Tuple{Int,Int}}, nmap::Dict{Int,Int})
    es = Tuple{Int,Int}[]
    for (u, v) in edges
        push!(es, (nmap[u], nmap[v]))
    end
    return es
end

function relabel(neighbormap::Vector{Dict{Int,Int}},
                 nmap::Dict{Int,Int}, emap::Dict{Int,Int})
    ns = Dict{Int,Int}[]
    for nbr in neighbormap
        d = Dict{Int,Int}()
        for (inc, adj) in nbr
            d[emap[inc]] = nmap[adj]
        end
        push!(ns, d)
    end
    return ns
end
