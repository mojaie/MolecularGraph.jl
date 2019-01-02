#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    VF2NodeInducedState,
    is_isomorphic,
    is_subgraph,
    isomorphmap,
    isomorphmapiter,
    vf2match!,
    updatestate!,
    candidatepairs,
    is_feasible,
    is_semantic_feasible,
    restore!


mutable struct VF2NodeInducedState <: VF2State
    G::UDGraph
    H::UDGraph

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

    function VF2NodeInducedState(G, H, mode, lim, nmatch, ematch, mand, forb)
        new(G, H, mode, lim, nmatch, ematch, mand, forb,
            Dict(), Dict(), Dict(), Dict())
    end
end


function is_isomorphic(G, H; kwargs...)
    return isomorphmap(G, H; mode=:graph, kwargs...) !== nothing
end


function is_subgraph(G, H; kwargs...)
    """ True if G is an induced subgraph of H"""
    return isomorphmap(H, G; kwargs...) !== nothing
end


function isomorphmap(G::UDGraph, H::UDGraph;
                     mode=:subgraph, depthlimit=1000,
                     nodematcher=nothing, edgematcher=nothing,
                     mandatory=nothing, forbidden=nothing)
    if nodecount(G) == 0 || nodecount(H) == 0
        return
    end
    state = VF2NodeInducedState(
        G, H, mode, depthlimit, nodematcher, edgematcher,
        mandatory, forbidden)
    return vf2match!(state)
end


function isomorphmapiter(G::UDGraph, H::UDGraph;
                         mode=:subgraph, depthlimit=1000,
                         nodematcher=nothing, edgematcher=nothing,
                         mandatory=nothing, forbidden=nothing)
    if nodecount(G) == 0 || nodecount(H) == 0
        return
    end
    state = VF2NodeInducedState(
        G, H, mode, depthlimit, nodematcher, edgematcher,
        mandatory, forbidden)
    return Channel(c::Channel -> vf2match!(state, c), ctype=Dict{Int,Int})
end


function vf2match!(state::VF2State, channel=nothing)
    # Recursive
    # println("depth $(length(state.g_core))")
    if length(state.g_core) == nodecount(state.H)
        # println("done $(state.g_core)")
        if channel === nothing
            return copy(state.g_core)
        else
            put!(channel, copy(state.g_core))
            return
        end
    end
    if length(state.g_core) >= state.depthlimit
        throw(OperationError("Maximum recursion reached"))
    end
    for (g, h) in candidatepairs(state)
        # println("candidates $(g) $(h)")
        if !is_feasible(state, g, h) || !is_semantic_feasible(state, g, h)
            continue
        end
        updatestate!(state, g, h)
        # println("g_core $(state.g_core)")
        mp = vf2match!(state, channel)
        if mp !== nothing
            return mp
        end
        # println("restored $(state.g_core)")
        restore!(state, g, h)
    end
    return
end


function updatestate!(state::VF2State, g, h)
    state.g_core[g] = h
    state.h_core[h] = g
    depth = length(state.g_core)
    if !haskey(state.g_term, g)
        state.g_term[g] = depth
    end
    if !haskey(state.h_term, h)
        state.h_term[h] = depth
    end
    g_nbrset = union([neighborkeys(state.G, n) for n in keys(state.g_core)]...)
    for n in setdiff(g_nbrset, keys(state.g_term))
        state.g_term[n] = depth
    end
    h_nbrset = union([neighborkeys(state.H, n) for n in keys(state.h_core)]...)
    for n in setdiff(h_nbrset, keys(state.h_term))
        state.h_term[n] = depth
    end
    return
end


function candidatepairs(state::VF2State)
    if state.mandatory !== nothing
        # Mandatory pair
        md = setdiff(keys(state.mandatory), keys(state.g_core))
        if !isempty(md)
            n = pop!(md)
            return [(n, state.mandatory[n])]
        end
    end

    pairs = Tuple{Int,Int}[]
    g_cand = setdiff(keys(state.g_term), keys(state.g_core))
    h_cand = setdiff(keys(state.h_term), keys(state.h_core))
    if isempty(g_cand) || isempty(h_cand)
        # New connected component
        g_cand = setdiff(nodekeys(state.G), keys(state.g_core))
        h_cand = setdiff(nodekeys(state.H), keys(state.h_core))
    end
    if !isempty(h_cand)
        h_min = minimum(h_cand)
        for g in g_cand
            if state.forbidden !== nothing
                # Forbidden pair
                if haskey(state.forbidden, g) && state.forbidden[g] == h_min
                    continue
                end
            end
            push!(pairs, (g, h_min))
        end
    end
    return pairs
end


function is_feasible(state::VF2State, g, h)
    # assume no self loop
    # Neighbor connectivity
    g_nbrs = neighborkeys(state.G, g)
    h_nbrs = neighborkeys(state.H, h)
    for n in intersect(g_nbrs, keys(state.g_core))
        if !(state.g_core[n] in h_nbrs)
            return false
        end
    end
    for n in intersect(h_nbrs, keys(state.h_core))
        if !(state.h_core[n] in g_nbrs)
            return false
        end
    end
    # Terminal set size
    g_term_count = length(setdiff(keys(state.g_term), keys(state.g_core)))
    h_term_count = length(setdiff(keys(state.h_term), keys(state.h_core)))
    if state.mode == :graph && g_term_count != h_term_count
        return false
    elseif state.mode == :subgraph && g_term_count < h_term_count
        return false
    end
    # Yet unexplored size
    g_new_count = length(setdiff(nodekeys(state.G), keys(state.g_term)))
    h_new_count = length(setdiff(nodekeys(state.H), keys(state.h_term)))
    if state.mode == :graph && g_new_count != h_new_count
        return false
    elseif state.mode == :subgraph && g_new_count < h_new_count
        return false
    end
    return true
end


function is_semantic_feasible(state::VF2NodeInducedState, g, h)
    if state.nodematch !== nothing
        if !state.nodematch(g, h)
            return false
        end
    end
    if state.edgematch !== nothing
        for nbr in intersect(neighborkeys(state.G, g), keys(state.g_core))
            g_edge = neighbors(state.G, g)[nbr]
            h_edge = neighbors(state.H, h)[state.g_core[nbr]]
            if !state.edgematch(g_edge, h_edge)
                return false
            end
        end
    end
    return true
end


function restore!(state::VF2State, g, h)
    depth = length(state.g_core)
    if g !== nothing && h !== nothing
        delete!(state.g_core, g)
        delete!(state.h_core, h)
    end
    for (k, v) in state.g_term
        if v == depth
            delete!(state.g_term, k)
        end
    end
    for (k, v) in state.h_term
        if v == depth
            delete!(state.h_term, k)
        end
    end
end
