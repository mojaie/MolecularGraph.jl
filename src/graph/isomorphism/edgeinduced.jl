#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#



function delta_y_correction!(mapping, G, H)
    gsub = edgesubgraph(G, Set(keys(mapping)))
    hsub = edgesubgraph(H, Set(values(mapping)))
    revmap = Dict(v => k for (k, v) in mapping)
    for e in delta_y_edges(mapping, gsub, hsub)
        delete!(mapping, e)
    end
    for e in delta_y_edges(revmap, hsub, gsub)
        delete!(mapping, revmap[e])
    end
end


function delta_y_edges(mapping, gsub, hsub)
    dys = Int[]
    for gn in triangles(gsub)
        g_edges = edgeset(nodesubgraph(gsub, Set(gn)))
        h_edges = Set([mapping[e] for e in g_edges])
        if nodecount(edgesubgraph(hsub, h_edges)) != 3
            push!(dys, pop!(g_edges))
        end
    end
    return dys
end


function delta_y_mismatch(G, H, mapping)
    # Delta-Y check for edge induced subgraph isomorphism
    gsub = edgesubgraph(G, Set(keys(mapping)))
    hsub = edgesubgraph(H, Set(values(mapping)))
    revmap = Dict(v => k for (k, v) in mapping)
    return delta_y(gsub, hsub, mapping) || delta_y(hsub, gsub, revmap)
end


function delta_y(gsub, hsub, mapping)
    for gn in triangles(gsub)
        g_edges = edgeset(nodesubgraph(gsub, Set(gn)))
        h_edges = Set([mapping[e] for e in g_edges])
        if nodecount(edgesubgraph(hsub, h_edges)) != 3
            return true
        end
    end
    return false
end


function lgnodematcher(G::LineGraph, H::LineGraph,
                       nodematcher::Function, edgematcher::Function)
    return function (g, h)
        edgematcher(g, h) || return false
        ge = nodeattr(G, g)
        he = nodeattr(H, h)
        m1 = nodematcher(ge.n1, he.n1) && nodematcher(ge.n2, he.n2)
        m2 = nodematcher(ge.n1, he.n2) && nodematcher(ge.n2, he.n1)
        return m1 || m2
    end
end


function lgedgematcher(G::LineGraph, H::LineGraph, nodematcher::Function)
    return (g, h) -> nodematcher(edgeattr(G, g).node, edgeattr(H, h).node)
end
