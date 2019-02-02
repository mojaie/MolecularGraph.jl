#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    isomorphismvf2,
    isomorphismitervf2,
    edgeisomorphismvf2,
    edgeisomorphismitervf2


mutable struct VF2State
    # Input
    G::UDGraph
    H::UDGraph
    mode::Symbol
    subgraphtype::Symbol
    # Optional
    timeout # Int
    nodematcher # Function
    edgematcher # Function
    mandatory # Dict{Int,Int}
    forbidden # Dict{Int,Int}
    # States
    g_core::Dict{Int,Int}
    h_core::Dict{Int,Int}
    g_term::Dict{Int,Int}
    h_term::Dict{Int,Int}
    expire # UInt64, nanoseconds
    status::Symbol

    function VF2State(G, H, mode, subg)
        state = new()
        state.G = G
        state.H = H
        state.mode = mode
        state.subgraphtype = subg
        state.g_core = Dict()
        state.h_core = Dict()
        state.g_term = Dict()
        state.h_term = Dict()
        state.status = :Ready
        return state
    end
end


function isomorphismvf2(G, H; kwargs...)
    return iterate(isomorphismitervf2(G, H; kwargs...))
end


function isomorphismitervf2(G, H; kwargs...)
    if nodecount(G) == 0 || nodecount(H) == 0
        return ()
    end
    state = VF2State(G, H, kwargs[:mode], kwargs[:subgraphtype])
    if haskey(kwargs, :timeout)
        state.timeout = kwargs[:timeout]::Int
        state.expire = (time_ns() + state.timeout * 1_000_000_000)::UInt64
    end
    (haskey(kwargs, :nodematcher)
        && (state.nodematcher = kwargs[:nodematcher]::Function))
    (haskey(kwargs, :edgematcher)
        && (state.edgematcher = kwargs[:edgematcher]::Function))
    (haskey(kwargs, :mandatory)
        && (state.mandatory = kwargs[:mandatory]::Dict{Int,Int}))
    (haskey(kwargs, :forbidden)
        && (state.forbidden = kwargs[:forbidden]::Dict{Int,Int}))
    return Channel(c::Channel -> expand!(state, c), ctype=Dict{Int,Int})
end


function edgeisomorphismvf2(G, H; kwargs...)
    return iterate(edgeisomorphismitervf2(G, H; kwargs...))
end


function edgeisomorphismitervf2(G, H;
        nodematcher=(g,h)->true, edgematcher=(g,h)->true, kwargs...)
    lg = linegraph(G)
    lh = linegraph(H)
    nmatch = lgnodematcher(lg, lh, nodematcher, edgematcher)
    ematch = lgedgematcher(lg, lh, nodematcher)
    results = isomorphismitervf2(lg, lh,
        nodematcher=nmatch, edgematcher=ematch; kwargs...)
    return Iterators.filter(results) do mapping
        return !delta_y_mismatch(G, H, mapping)
    end
end


function expand!(state::VF2State, channel)
    # Recursive
    # println("depth $(length(state.g_core))")
    if length(state.g_core) == nodecount(state.H)
        # println("done $(state.g_core)")
        put!(channel, copy(state.g_core))
        return
    elseif isdefined(state, :timeout) && time_ns() > state.expire
        state.status = :Timedout
        return
    end
    for (g, h) in candidatepairs(state)
        # println("candidates $(g) $(h)")
        if !is_feasible(state, g, h) || !is_semantic_feasible(state, g, h)
            continue
        end
        updatestate!(state, g, h)
        # println("g_core $(state.g_core)")
        mp = expand!(state, channel)
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


function candidatepairs(state::VF2State)
    if isdefined(state, :mandatory)
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
            if isdefined(state, :forbidden)
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
    if state.mode == :Isomorphism && g_term_count != h_term_count
        return false
    elseif state.mode == :Subgraph && g_term_count < h_term_count
        return false
    end
    # Yet unexplored size
    g_new_count = length(setdiff(nodekeys(state.G), keys(state.g_term)))
    h_new_count = length(setdiff(nodekeys(state.H), keys(state.h_term)))
    if state.mode == :Isomorphism && g_new_count != h_new_count
        return false
    elseif state.mode == :Subgraph && g_new_count < h_new_count
        return false
    end
    return true
end


function is_semantic_feasible(state::VF2State, g, h)
    if isdefined(state, :nodematcher)
        if !state.nodematcher(g, h)
            return false
        end
    end
    if isdefined(state, :edgematcher)
        for nbr in intersect(neighborkeys(state.G, g), keys(state.g_core))
            g_edge = neighbors(state.G, g)[nbr]
            h_edge = neighbors(state.H, h)[state.g_core[nbr]]
            if !state.edgematcher(g_edge, h_edge)
                return false
            end
        end
    end
    return true
end
