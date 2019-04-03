#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    LineGraph, linegraph


struct LineGraphNode <: AbstractNode
    n1::Int
    n2::Int
end


struct LineGraphEdge <: UndirectedEdge
    node::Int
end


struct LineGraph <: OrderedGraph
    neighbormap::Vector{Dict{Int,Int}}
    edges::Vector{Tuple{Int,Int}}
    nodeattrs::Vector{LineGraphNode}
    edgeattrs::Vector{LineGraphEdge}
    cache::Dict{Symbol,Any}
end


"""
    linegraph(G::AbstractGraph) -> LineGraph

Generate line graph.
"""
function linegraph(G::AbstractGraph)
    L = LineGraph([], [], [], [], Dict())
    nmap = Dict{Int,Int}()
    for (i, e) in enumerate(edgeset(G))
        (u, v) = getedge(G, e)
        push!(L.nodeattrs, LineGraphNode(u, v))
        push!(L.neighbormap, Dict())
        nmap[e] = i
    end
    ecnt = 1
    for n in nodeset(G)
        degree(G, n) < 2 && continue
        for (e1, e2) in combinations(incidences(G, n))
            ne1 = nmap[e1]
            ne2 = nmap[e2]
            push!(L.edges, (ne1, ne2))
            push!(L.edgeattrs, LineGraphEdge(n))
            L.neighbormap[ne1][ecnt] = ne2
            L.neighbormap[ne2][ecnt] = ne1
            ecnt += 1
        end
    end
    return L
end

function linegraph(G::OrderedGraph)
    L = LineGraph([], [], [], [], Dict())
    for (i, edge) in enumerate(edgesiter(G))
        push!(L.nodeattrs, LineGraphNode(edge...))
        push!(L.neighbormap, Dict())
    end
    ecnt = 1
    for i in 1:nodecount(G)
        degree(G, i) < 2 && continue
        for (e1, e2) in combinations(incidences(G, i))
            push!(L.edges, (e1, e2))
            push!(L.edgeattrs, LineGraphEdge(i))
            L.neighbormap[e1][ecnt] = e2
            L.neighbormap[e2][ecnt] = e1
            ecnt += 1
        end
    end
    return L
end
