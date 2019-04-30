#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    ModularProduct, modularproduct,
    CartesianProduct, cartesianproduct


# Modular product

struct ModularProductNode <: AbstractNode
    g::Int
    h::Int
end


struct ModularProductEdge <: UndirectedEdge
    hasedge::Bool
end


struct ModularProduct <: OrderedGraph
    neighbormap::Vector{Dict{Int,Int}}
    edges::Vector{Tuple{Int,Int}}
    nodeattrs::Vector{ModularProductNode}
    edgeattrs::Vector{ModularProductEdge}
    cache::Dict{Symbol,Any}
end



"""
    modularproduct(G::OrderedGraph, H::OrderedGraph) -> ModularProduct

Return the modular product of graphs G and H.
"""
function modularproduct(G::OrderedGraph, H::OrderedGraph;
            nodematcher=(g,h)->true,
            edgefilter=(g1,g2,h1,h2)->hasedge(G,g1,g2)==hasedge(H,h1,h2),
            kwargs...)
    product = ModularProduct([], [], [], [], Dict())
    ndict = Dict{Int,Dict{Int,Int}}() # Ref to node indices of the product
    # Modular product nodes
    nattrs = ModularProductNode[]
    for g in 1:nodecount(G)
        ndict[g] = Dict{Int,Int}()
        for h in 1:nodecount(H)
            ndict[g][h] = addnode!(product, ModularProductNode(g, h))
        end
    end
    (nodecount(G) < 2 || nodecount(H) < 2) && return product
    # Modular product edges
    for (g1, g2) in combinations(1:nodecount(G))
        for (h1, h2) in combinations(1:nodecount(H))
            edgefilter(g1, g2, h1, h2) || continue
            if nodematcher(g1, h1) && nodematcher(g2, h2)
                addedge!(
                    product, ndict[g1][h1], ndict[g2][h2],
                    ModularProductEdge(hasedge(G, g1, g2))
                )
            end
            if nodematcher(g1, h2) && nodematcher(g2, h1)
                addedge!(
                    product, ndict[g1][h2], ndict[g2][h1],
                    ModularProductEdge(hasedge(G, g1, g2))
                )
            end
        end
    end
    return product
end



# Cartesian product

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
