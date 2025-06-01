#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    twocoloring, triangles, cartesianproduct


"""
    twocoloring(graph::UndirectedGraph)

Do 2-coloring and return two sets of nodes.
"""
function twocoloring(graph::UndirectedGraph)
    nodes = nodeset(graph)
    c1 = Set()
    c2 = Set()
    while !isempty(nodes)
        root = pop!(nodes)
        queue = [root]
        push!(c1, root)
        while !isempty(queue)
            i = popfirst!(queue)
            for adj in adjacencies(graph, i)
                if !(adj in c1 || adj in c2)
                    if i in c1
                        push!(c2, adj)
                        push!(queue, adj)
                    else
                        push!(c1, adj)
                        push!(queue, adj)
                    end
                elseif i in c1 && adj in c1
                    return  # No 2-coloring
                elseif i in c2 && adj in c2
                    return  # No 2-coloring
                end
            end
        end
        setdiff!(nodes, c1, c2)
    end
    return c1, c2
end


"""
    triangles(graph::UndirectedGraph) -> Set{Tuple{Int,Int,Int}}

Fast computation of finding triangle node sets in the graph.
"""
function triangles(graph::UndirectedGraph)
    # deprecated: better implementation in Graphs.localclustering!
    triads = Set()
    for n in nodeset(graph)
        adjs = collect(adjacencies(graph, n))
        length(adjs) < 2 && continue
        for (e1, e2) in combinations(length(adjs))
            u, v = (adjs[e1], adjs[e2])
            if hasedge(graph, u, v)
                push!(triads, tuple(sort([n, u, v])...))
            end
        end
    end
    return triads
end


# Cartesian product
#
struct CartesianProductNode <: AbstractNode
    g::Int
    h::Int
end


struct CartesianProduct <: OrderedGraph
    neighbormap::Vector{Dict{Int,Int}}
    edges::Vector{Tuple{Int,Int}}
    nodeattrs::Vector{CartesianProductNode}
    cache::Dict{Symbol,Any}
end


"""
    cartesianproduct(G::OrderedGraph, H::OrderedGraph) -> CartesianProduct

Return the cartesian product of graphs G and H.
"""
function cartesianproduct(G::OrderedGraph, H::OrderedGraph)
    product = CartesianProduct([], [], [], Dict())
    ndict = Dict{Int,Dict{Int,Int}}() # Ref to node indices of the product
    #Cartesian product nodes
    nattrs = CartesianProductNode[]
    for g in 1:nodecount(G)
        ndict[g] = Dict{Int,Int}()
        for h in 1:nodecount(H)
            ndict[g][h] = addnode!(product, CartesianProductNode(g, h))
        end
    end
    (nodecount(G) < 2 || nodecount(H) < 2) && return product
    # Cartesian product edges
    for i in 1:nodecount(G)
        for (u, v) in edgesiter(H)
            addedge!(product, ndict[i][u], ndict[i][v])
        end
    end
    for i in 1:nodecount(H)
        for (u, v) in edgesiter(G)
            addedge!(product, ndict[u][i], ndict[v][i])
        end
    end
    return product
end
