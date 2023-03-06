#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    ModularProduct, modular_product


# Modular product

struct ModularProduct{T}
    graph::SimpleGraph{T}
    isconnected::Dict{Edge{T},Bool}
end

# modular_product_reverse_map
# i -> (div(i, nv(h)), mod(i, nv(h)))

"""
    modularproduct(G::OrderedGraph, H::OrderedGraph) -> ModularProduct

Return the modular product of graphs G and H.
"""
function modular_product(g::SimpleGraph{T}, h::SimpleGraph{T};
            nodematcher=(g1,h1)->true,
            edgefilter=(g1,g2,h1,h2)->has_edge(g,g1,g2)==has_edge(g,h1,h2),
            kwargs...) where T
    m = SimpleGraph(nv(g) * nv(h))
    connected = Dict{Edge{T},Bool}()
    id(i, j) = (i - 1) * nv(h) + j
    for (g1, g2) in combinations(nv(g))
        for (h1, h2) in combinations(nv(h))
            edgefilter(g1, g2, h1, h2) || continue
            if nodematcher(g1, h1) && nodematcher(g2, h2)
                e = Edge{T}(id(g1, h1), id(g2, h2))  # g1 < g2 -> id(g1, h1) < id(g2, h2)
                add_edge!(m, e)  
                connected[e] = has_edge(g, g1, g2)
            end
            if nodematcher(g1, h2) && nodematcher(g2, h1)
                e = Edge{T}(id(g1, h2), id(g2, h1))  # g1 < g2 -> id(g1, h2) < id(g2, h1)
                add_edge!(m, e)
                connected[e] = has_edge(g, g1, g2)
            end
        end
    end
    return ModularProduct(m, connected)
end
