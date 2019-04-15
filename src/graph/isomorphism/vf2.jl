#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    isomorphismitervf2, edgeisomorphismitervf2,
    graphmatches, graphmatch, isgraphmatch,
    subgraphmatches, subgraphmatch, issubgraphmatch,
    edgesubgraphmatches, edgesubgraphmatch, isedgesubgraphmatch


mutable struct VF2State{T1<:UndirectedGraph,T2<:UndirectedGraph}
    G::T1
    H::T2
    mode::Symbol
    timeout # Int
    nodematcher # Function
    edgematcher # Function
    mandatory # Dict{Int,Int}
    forbidden # Dict{Int,Int}

    expire # UInt64, nanoseconds

    g_core::Dict{Int,Int}
    h_core::Dict{Int,Int}
    g_term::Dict{Int,Int}
    h_term::Dict{Int,Int}
    status::Symbol

    function VF2State{T1,T2}(G::T1, H::T2, mode
            ) where {T1<:UndirectedGraph,T2<:UndirectedGraph}
        state = new()
        state.G = G
        state.H = H
        state.mode = mode
        state.g_core = Dict()
        state.h_core = Dict()
        state.g_term = Dict()
        state.h_term = Dict()
        state.status = :Ready
        return state
    end
end


function expand!(state::VF2State, channel; verbose=false, kwargs...)
    # Recursive
    verbose && println("depth $(length(state.g_core))")
    if length(state.g_core) == nodecount(state.H)
        verbose && println("done $(state.g_core)")
        put!(channel, copy(state.g_core))
        return
    elseif isdefined(state, :timeout) && time_ns() > state.expire
        state.status = :Timedout
        verbose && println("Timed out")
        return
    end
    for (g, h) in candidatepairs(state)
        verbose && println("candidates $(g) $(h)")
        if (!is_feasible(state, g, h, verbose=verbose)
                || !is_semantic_feasible(state, g, h, verbose=verbose))
            continue
        end
        updatestate!(state, g, h)
        verbose && println("g_core $(state.g_core)")
        mp = expand!(state, channel, verbose=verbose)
        verbose && println("restored $(state.g_core)")
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
    g_nbrset = union([adjacencies(state.G, n) for n in keys(state.g_core)]...)
    for n in setdiff(g_nbrset, keys(state.g_term))
        state.g_term[n] = depth
    end
    h_nbrset = union([adjacencies(state.H, n) for n in keys(state.h_core)]...)
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
        g_cand = setdiff(nodeset(state.G), keys(state.g_core))
        h_cand = setdiff(nodeset(state.H), keys(state.h_core))
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


function is_feasible(state::VF2State, g, h; verbose=false)
    # TODO: assume no self loop
    g_nbrs = adjacencies(state.G, g)
    h_nbrs = adjacencies(state.H, h)
    for n in intersect(g_nbrs, keys(state.g_core))
        if !(state.g_core[n] in h_nbrs)
            verbose && println("Infeasible: neighbor connectivity")
            return false
        end
    end
    for n in intersect(h_nbrs, keys(state.h_core))
        if !(state.h_core[n] in g_nbrs)
            verbose && println("Infeasible: neighbor connectivity")
            return false
        end
    end
    g_term_count = length(setdiff(keys(state.g_term), keys(state.g_core)))
    h_term_count = length(setdiff(keys(state.h_term), keys(state.h_core)))
    if state.mode == :Isomorphism && g_term_count != h_term_count
        verbose && println("Infeasible: terminal set size")
        return false
    elseif state.mode == :Subgraph && g_term_count < h_term_count
        verbose && println("Infeasible: terminal set size")
        return false
    end
    g_new_count = length(setdiff(nodeset(state.G), keys(state.g_term)))
    h_new_count = length(setdiff(nodeset(state.H), keys(state.h_term)))
    if state.mode == :Isomorphism && g_new_count != h_new_count
        verbose && println("Infeasible: yet unexplored nodes")
        return false
    elseif state.mode == :Subgraph && g_new_count < h_new_count
        verbose && println("Infeasible: yet unexplored nodes")
        return false
    end
    verbose && println("Syntactically feasible")
    return true
end


function is_semantic_feasible(state::VF2State, g, h; verbose=false)
    if isdefined(state, :nodematcher)
        if !state.nodematcher(g, h)
            verbose && println("Infeasible: node attribute mismatch")
            return false
        end
    end
    if isdefined(state, :edgematcher)
        for (inc, adj) in neighbors(state.G, g)
            haskey(state.g_core, adj) || continue
            hinc = findedgekey(state.H, h, state.g_core[adj])
            if !state.edgematcher(inc, hinc)
                verbose && println("Infeasible: edge attribute mismatch")
                return false
            end
        end
    end
    return true
end



function isomorphismitervf2(G::T1, H::T2; kwargs...
        ) where {T1<:UndirectedGraph,T2<:UndirectedGraph}
    (nodecount(G) == 0 || nodecount(H) == 0) && return ()
    state = VF2State{T1,T2}(G, H, kwargs[:mode])
    if haskey(kwargs, :nodematcher)
        state.nodematcher = kwargs[:nodematcher]::Function
    end
    if haskey(kwargs, :edgematcher)
        state.edgematcher = kwargs[:edgematcher]::Function
    end
    if haskey(kwargs, :mandatory)
        state.mandatory = kwargs[:mandatory]::Dict{Int,Int}
    end
    if haskey(kwargs, :forbidden)
        state.forbidden = kwargs[:forbidden]::Dict{Int,Int}
    end
    if haskey(kwargs, :timeout)
        state.timeout = kwargs[:timeout]::Int
        state.expire = (time_ns() + state.timeout * 1_000_000_000)::UInt64
    end
    return Channel(
        c::Channel -> expand!(state, c; kwargs...), ctype=Dict{Int,Int})
end



function edgeisomorphismitervf2(G::T1, H::T2; kwargs...
        ) where {T1<:UndirectedGraph,T2<:UndirectedGraph}
    (edgecount(G) == 0 || edgecount(H) == 0) && return ()
    lg = linegraph(G)
    lh = linegraph(H)
    state = VF2State{LineGraph,LineGraph}(lg, lh, kwargs[:mode])
    if haskey(kwargs, :nodematcher)
        state.edgematcher = lgedgematcher(lg, lh, kwargs[:nodematcher])
        if haskey(kwargs, :edgematcher)
            state.nodematcher = lgnodematcher(
                lg, lh, kwargs[:nodematcher], kwargs[:edgematcher])
        end
    end
    if haskey(kwargs, :mandatory)
        state.mandatory = kwargs[:mandatory]::Dict{Int,Int}
    end
    if haskey(kwargs, :forbidden)
        state.forbidden = kwargs[:forbidden]::Dict{Int,Int}
    end
    if haskey(kwargs, :timeout)
        state.timeout = kwargs[:timeout]::Int
        state.expire = (time_ns() + state.timeout * 1_000_000_000)::UInt64
    end
    results = Channel(
        c::Channel -> expand!(state, c; kwargs...), ctype=Dict{Int,Int})
    return Iterators.filter(results) do mapping
        return !delta_y_mismatch(G, H, mapping)
    end
end



# Graph isomorphism

"""
    graphmatches(G::AbstractGraph, H::AbstractGraph; kwargs...) -> Iterator

Generate isomorphism mappings between `G` and `H`. If no match found, return
nothing.
"""
graphmatches(G, H; kwargs...
    ) = isomorphismitervf2(G, H, mode=:Isomorphism; kwargs...)

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
    subgraphmatches(G::AbstractGraph, H::AbstractGraph; kwargs...) -> Iterator

Generate subgraph isomorphism mappings between `G` and `H`.

# Keyword arguments

- nodematcher(Function): node matcher function that takes two node indices as
arguments.
- edgematcher(Function): edge matcher function that takes two edge indices as
arguments.
- mandatory(Dict{Int,Int}): mandatory node matches (available for only VF2)
- forbidden(Dict{Int,Int}):
    forbidden node matches (available for only VF2)
"""
subgraphmatches(G, H; kwargs...
    ) = isomorphismitervf2(G, H, mode=:Subgraph; kwargs...)

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
edgesubgraphmatches(G, H; kwargs...
    ) = edgeisomorphismitervf2(G, H, mode=:Subgraph; kwargs...)

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
