#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    ModularProduct,
    modularproduct,
    edgefilter


struct ModularProductNode <: AbstractNode
    g::Int
    h::Int
end


struct ModularProductEdge <: UndirectedEdge
    u::Int
    v::Int
    hasedge::Bool
end


ModularProduct = MapUDGraph{ModularProductNode,ModularProductEdge}


"""
    modularproduct(G::UDGraph, H::UDGraph) -> MapUDGraph

Return the modular product of graphs G and H.
"""
function modularproduct(G::UDGraph, H::UDGraph,
        nodematcher=(g,h)->true, edgefilter=edgefilter(G, H))
    product = ModularProduct()
    ncnt = 0
    ndict = Dict{Int,Dict{Int,Int}}() # Ref to node indices of the product
    # Modular product nodes
    for g in nodekeys(G)
        ndict[g] = Dict{Int,Int}()
        for h in nodekeys(H)
            ncnt += 1
            n = ModularProductNode(g, h)
            updatenode!(product, n, ncnt)
            ndict[g][h] = ncnt
        end
    end
    if nodecount(G) < 2 || nodecount(H) < 2
        return product
    end
    # Modular product edges
    ecnt = 0
    for (g1, g2) in combinations(nodekeys(G))
        for (h1, h2) in combinations(nodekeys(H))
            if !edgefilter(g1, g2, h1, h2)
                continue
            end
            if nodematcher(g1, h1) && nodematcher(g2, h2)
                # TODO: hasedge
                e = ModularProductEdge(
                    ndict[g1][h1], ndict[g2][h2], g2 in neighborkeys(G, g1))
                ecnt += 1
                updateedge!(product, e, ecnt)
            end
            if nodematcher(g1, h2) && nodematcher(g2, h1)
                e = ModularProductEdge(
                    ndict[g1][h2], ndict[g2][h1], g2 in neighborkeys(G, g1))
                ecnt += 1
                updateedge!(product, e, ecnt)
            end
        end
    end
    return product
end


function edgefilter(G, H)
    return function (g1, g2, h1, h2)
        return (g2 in neighborkeys(G, g1)) == (h2 in neighborkeys(H, h1))
    end
end
