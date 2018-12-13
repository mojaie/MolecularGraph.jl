#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    VF2EdgeInducedState,
    vf2edgesubgraphstate,
    edge_subgraph_isomorph,
    isomorphmap!,
    is_semantic_feasible,
    deltaYmismatch


mutable struct VF2EdgeInducedState <: VF2State
    G::AbstractUGraph
    H::AbstractUGraph

    mode::Symbol
    depthlimit::Int
    nodematch::Union{Function,Nothing}
    edgematch::Union{Function,Nothing}
    mandatory::Dict{Int,Int}
    forbidden::Dict{Int,Int}

    g_core::Dict{Int,Int}
    h_core::Dict{Int,Int}
    g_term::Dict{Int,Int}
    h_term::Dict{Int,Int}

    mappings::Vector{Dict{Int,Int}}
end

vf2edgesubgraphstate(G, H) = VF2EdgeInducedState(
    G, H, :subgraph, 1000, nothing, nothing, Dict(), Dict(),
    Dict(), Dict(), Dict(), Dict(), []
)


function edge_subgraph_isomorph(G, H)
    """ True if H is an edge-induced subgraph of G"""
    state = vf2edgesubgraphstate(G, H)
    isomorphmap!(state)
    !isempty(state.mappings)
end


function isomorphmap!(state::VF2EdgeInducedState)
    # Edge induced subgraph isomorphism mapping
    G = state.G
    H = state.H
    state.G = linegraph(G)
    state.H = linegraph(H)
    vf2match!(state, nothing, nothing)
    deltaYfiltered = []
    for mapping in state.mappings
        if !deltaYmismatch(G, H, mapping)
            push!(deltaYfiltered, mapping)
        end
    end
    state.mappings = deltaYfiltered
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


function deltaYmismatch(G, H, mapping)
    # Delta-Y check for edge induced subgraph isomorphism
    gsub = edgesubgraph(G, keys(mapping))
    hsub = edgesubgraph(H, values(mapping))
    revmap = Dict(v => k for (k, v) in mapping)
    return deltaY(gsub, hsub, mapping) || deltaY(hsub, gsub, revmap)
end


function deltaY(gsub, hsub, mapping)
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
