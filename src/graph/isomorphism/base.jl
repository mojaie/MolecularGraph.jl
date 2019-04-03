#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    graphisomorphism,
    graphmatches, graphmatch, isgraphmatch,
    subgraphmatches, subgraphmatch, issubgraphmatch,
    edgesubgraphmatches, edgesubgraphmatch, isedgesubgraphmatch,
    mcismap, mcissize, mcesmap, mcessize


"""
    graphisomorphism(G::AbstractGraph, H::AbstractGraph;
        mode=:Subgraph, algorithm=:VF2, subgraphtype=:NodeInduced,
        kwargs...) -> Bool

Compute (sub)graph isomorphism.

## Required parameters:

- algorithm: `:VF2`, `:Clique`
- mode: `:Isomorphism`(VF2), `:Subgraph`(VF2), `:MCS`(Clique)
- subgraphtype: `:NodeInduced`, `:EdgeInduced` (for mode=:Subgraph)
- connectivity: `:Connected`, `:Disconnected`  (for mode=:MCS)
- constraint: `TopologicalConstraint`, `DiameterRestriction` (for mode=:MCS)

## Optional parameters:

- timeout(Int):
    return suboptimal result when timed out (second).
- threshold(Int):
    return suboptimal result when the given threshold achieved.
- nodematcher(Function):
    node matcher function that takes two node indices as arguments.
- edgematcher(Function):
    edge matcher function that takes two edge indices as arguments.
- mandatory(Dict{Int,Int}):
    mandatory node matches (available for only VF2)
- forbidden(Dict{Int,Int}):
    forbidden node matches (available for only VF2)
- theta(Int):
    distance mismatch tolerance in topologically constrainted MCS
- diameter(Int):
    diameter size in MCS with diameter restriction
"""
function graphisomorphism(G::AbstractGraph, H::AbstractGraph;
        algorithm=:VF2, subgraphtype=:NodeInduced, mode=:Isomorphism, kwargs...)
    if algorithm == :VF2
        if subgraphtype == :NodeInduced
            return isomorphismitervf2(
                G, H, subgraphtype=:NodeInduced, mode=mode; kwargs...)
        elseif subgraphtype == :EdgeInduced
            return edgeisomorphismitervf2(
                G, H, subgraphtype=:EdgeInduced, mode=mode; kwargs...)
        end
        # TODO: implement MCS
    elseif algorithm == :Clique
        # TODO: implement subgraphmatch
        if mode == :MCS
            if subgraphtype == :NodeInduced
                return nodemcsclique(G, H, mode=mode; kwargs...)
            elseif subgraphtype == :EdgeInduced
                return edgemcsclique(G, H, mode=mode; kwargs...)
            end
        end
    end
end



# Graph isomorphism

"""
    graphmatches(G::AbstractGraph, H::AbstractGraph; kwargs...) -> Iterator

Generate isomorphism mappings between `G` and `H`. If no match found, return
nothing.
"""
graphmatches(G, H; kwargs...) = graphisomorphism(G, H; kwargs...)

"""
    graphmatch(G::AbstractGraph, H::AbstractGraph; kwargs...) -> Dict{Int,Int}

Return an isomorphism mapping between `G` and `H`. If no match found, return
nothing.
"""
function graphmatch(G, H; kwargs...)
    res = iterate(graphmatches(G, H; kwargs...))
    return res === nothing ? nothing : res[1]
end

"""
    isgraphmatch(G::AbstractGraph, H::AbstractGraph; kwargs...) -> Bool

Return true if `G` and `H` are isomorphic.
"""
isgraphmatch(G, H; kwargs...) = graphmatch(G, H; kwargs...) !== nothing



# Node induced subgraph isomorphism

"""
    subgraphmatch(G::AbstractGraph, H::AbstractGraph; kwargs...) -> Iterator

Generate subgraph isomorphism mappings between `G` and `H`.
"""
subgraphmatches(G, H; kwargs...
    ) = graphisomorphism(G, H, mode=:Subgraph; kwargs...)

"""
    subgraphmatch(
        G::AbstractGraph, H::AbstractGraph; kwargs...) -> Dict{Int,Int}

Return a subgraph isomorphism mapping between `G` and `H`. If no match found,
return nothing.
"""
function subgraphmatch(G, H; kwargs...)
    res = iterate(subgraphmatches(G, H; kwargs...))
    return res === nothing ? nothing : res[1]
end

"""
    issubgraphmatch(G::AbstractGraph, H::AbstractGraph; kwargs...) -> Bool

Return true if a node induced subgraph of `G` and `H` are isomorphic.
"""
issubgraphmatch(G, H; kwargs...) = subgraphmatch(G, H; kwargs...) !== nothing



# Edge induced subgraph isomorphism

"""
    edgesubgraphmatch(G::AbstractGraph, H::AbstractGraph; kwargs...) -> Iterator

Generate edge induced subgraph isomorphism mappings between `G` and `H`.
"""
edgesubgraphmatches(G, H; kwargs...) = graphisomorphism(
    G, H, mode=:Subgraph, subgraphtype=:EdgeInduced; kwargs...)

"""
    edgesubgraphmatch(
        G::AbstractGraph, H::AbstractGraph; kwargs...) -> Dict{Int,Int}

Return a edge induced subgraph isomorphism mapping between `G` and `H`.
If no match found, return nothing.
"""
function edgesubgraphmatch(G, H; kwargs...)
    res = iterate(edgesubgraphmatches(G, H; kwargs...))
    return res === nothing ? nothing : res[1]
end

"""
    isedgesubgraphmatch(G::AbstractGraph, H::AbstractGraph; kwargs...) -> Bool

Return true if a node induced subgraph of `G` and `H` are isomorphic.
"""
isedgesubgraphmatch(G, H; kwargs...
    ) = edgesubgraphmatch(G, H; kwargs...) !== nothing



# MCIS

"""
    mcismap(G::AbstractGraph, H::AbstractGraph; kwargs...) -> Dict{Int,Int}

Return a maximum common induced subgraph mapping of `G` and ``H.
"""
function mcismap(G, H; kwargs...)
    res = iterate(
        graphisomorphism(G, H, algorithm=:Clique, mode=:MCS; kwargs...))
    return res === nothing ? nothing : res[1]
end

"""
    mcissize(G::AbstractGraph, H::AbstractGraph; kwargs...) -> Int

Return the maximum common induced subgraph size (number of nodes).
"""
mcissize(G, H; kwargs...) = length(mcismap(G, H; kwargs...))


# MCES

"""
    mcesmap(G::AbstractGraph, H::AbstractGraph; kwargs...) -> Dict{Int,Int}

Return a maximum common edge induced subgraph mapping of `G` and ``H.
"""
function mcesmap(G, H; kwargs...)
    res = iterate(graphisomorphism(G, H,
        algorithm=:Clique, mode=:MCS, subgraphtype=:EdgeInduced; kwargs...))
    return res === nothing ? nothing : res[1]
end

"""
    mcessize(G::AbstractGraph, H::AbstractGraph; kwargs...) -> Int

Return the maximum common induced subgraph size (number of nodes).
"""
mcessize(G, H; kwargs...) = length(mcesmap(G, H; kwargs...))




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
