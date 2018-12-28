#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    VF2EdgeInducedState,
    is_edge_subgraph,
    edgeisomorph,
    edgeisomorphiter,
    is_semantic_feasible,
    delta_y_mismatch


mutable struct VF2EdgeInducedState <: VF2State
    G::AbstractUGraph
    H::AbstractUGraph

    mode::Symbol
    depthlimit::Int
    nodematch::Union{Function,Nothing}
    edgematch::Union{Function,Nothing}
    mandatory::Union{Dict{Int,Int},Nothing}
    forbidden::Union{Dict{Int,Int},Nothing}

    g_core::Dict{Int,Int}
    h_core::Dict{Int,Int}
    g_term::Dict{Int,Int}
    h_term::Dict{Int,Int}

    function VF2EdgeInducedState(G, H, mode, lim, nmatch, ematch, mand, forb)
        new(G, H, mode, lim, nmatch, ematch, mand, forb,
            Dict(), Dict(), Dict(), Dict())
    end
end


function is_edge_subgraph(G, H; kwargs...)
    """ True if G is an edge-induced subgraph of H"""
    return edgeisomorph(H, G; kwargs...) !== nothing
end


function edgeisomorph(G::AbstractUGraph, H::AbstractUGraph; kwargs...)
    return iterate(edgeisomorphiter(G, H; kwargs...))
end


function edgeisomorphiter(G::AbstractUGraph, H::AbstractUGraph;
                          mode=:subgraph, depthlimit=1000,
                          nodematcher=nothing, edgematcher=nothing,
                          mandatory=nothing, forbidden=nothing)
    # Edge induced subgraph isomorphism mapping
    if edgecount(G) == 0 || edgecount(H) == 0
        return ()
    end
    state = VF2EdgeInducedState(
        linegraph(G), linegraph(H), mode, depthlimit,
        nodematcher, edgematcher, mandatory, forbidden)
    channel = Channel(c::Channel -> vf2match!(state, c), ctype=Dict{Int,Int})
    return Iterators.filter(channel) do mapping
        return !delta_y_mismatch(G, H, mapping)
    end
end


function is_semantic_feasible(state::VF2EdgeInducedState, g, h)
    # TODO: refactor
    # Linegraph node match (edge match)
    if state.edgematch !== nothing
        if !state.edgematch(g, h)
            return false
        end
        gn = getnode(state.G, g)
        hn = getnode(state.H, h)
        m1 = state.nodematch(gn.n1, hn.n1) && state.nodematch(gn.n2, hn.n2)
        m2 = state.nodematch(gn.n2, hn.n1) && state.nodematch(gn.n1, hn.n2)
        if !m1 && !m2
            return false
        end
    end
    # Linegraph edge match (node match)
    if state.nodematch !== nothing
        for nbr in intersect(neighborkeys(state.G, g), keys(state.g_core))
            g_edge = neighbors(state.G, g)[nbr]
            h_edge = neighbors(state.H, h)[state.g_core[nbr]]
            gn = getedge(state.G, g_edge).node
            hn = getedge(state.H, h_edge).node
            if !state.nodematch(gn, hn)
                return false
            end
        end
    end
    return true
end


function delta_y_mismatch(G, H, mapping)
    # Delta-Y check for edge induced subgraph isomorphism
    gsub = edgesubgraph(G, keys(mapping))
    hsub = edgesubgraph(H, values(mapping))
    revmap = Dict(v => k for (k, v) in mapping)
    return delta_y(gsub, hsub, mapping) || delta_y(hsub, gsub, revmap)
end


function delta_y(gsub, hsub, mapping)
    for gn in triangles(gsub)
        triad = [(1, 2), (2, 3), (1, 3)]
        g_edges = [neighbors(gsub, gn[u])[gn[v]] for (u, v) in triad]
        h_edges = [mapping[e] for e in g_edges]
        if nodecount(edgesubgraph(hsub, h_edges)) != 3
            return true
        end
    end
    return false
end
