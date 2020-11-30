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

    g_core::Dict{Int,Int}
    h_core::Dict{Int,Int}
    g_term::Dict{Int,Int}
    h_term::Dict{Int,Int}
    expire::Union{UInt64,Nothing} # UInt64, nanoseconds

    mappings::Vector{Dict{Int,Int}}
    status::Symbol

    function VF2State{T1,T2}(G::T1, H::T2; timeout=nothing
            ) where {T1<:UndirectedGraph,T2<:UndirectedGraph}
        if timeout !== nothing
            expire = (time_ns() + timeout * 1_000_000_000)::UInt64
        else
            expire = nothing
        end
        return new(G, H, Dict(), Dict(), Dict(), Dict(), expire, [], :ready)
    end
end


function candidatepairs(state::VF2State; kwargs...)
    # Mandatory pair
    if haskey(kwargs, :mandatory)
        md = setdiff(keys(kwargs[:mandatory]), keys(state.g_core))
        if !isempty(md)
            n = pop!(md)
            return [(n, kwargs[:mandatory][n])]
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
            # Forbidden pair
            if haskey(kwargs, :forbidden)
                fb = kwargs[:forbidden]
                haskey(fb, g) && fb[g] == h_min && continue
            end
            push!(pairs, (g, h_min))
        end
    end
    return pairs
end


function is_feasible(state::VF2State, g, h; mode=:Subgraph, kwargs...)
    # TODO: assume no self loop
    g_nbrs = adjacencies(state.G, g)
    h_nbrs = adjacencies(state.H, h)
    for n in intersect(g_nbrs, keys(state.g_core))
        if !(state.g_core[n] in h_nbrs)
            # println("Infeasible: neighbor connectivity")
            return false
        end
    end
    for n in intersect(h_nbrs, keys(state.h_core))
        if !(state.h_core[n] in g_nbrs)
            # println("Infeasible: neighbor connectivity")
            return false
        end
    end
    g_term_count = length(setdiff(keys(state.g_term), keys(state.g_core)))
    h_term_count = length(setdiff(keys(state.h_term), keys(state.h_core)))
    if mode == :Isomorphism && g_term_count != h_term_count
        # println("Infeasible: terminal set size")
        return false
    elseif mode == :Subgraph && g_term_count < h_term_count
        # println("Infeasible: terminal set size")
        return false
    end
    g_new_count = length(setdiff(nodeset(state.G), keys(state.g_term)))
    h_new_count = length(setdiff(nodeset(state.H), keys(state.h_term)))
    if mode == :Isomorphism && g_new_count != h_new_count
        # println("Infeasible: yet unexplored nodes")
        return false
    elseif mode == :Subgraph && g_new_count < h_new_count
        # println("Infeasible: yet unexplored nodes")
        return false
    end
    # println("Syntactically feasible")
    return true
end


function is_semantic_feasible(state::VF2State, g, h; kwargs...)
    if haskey(kwargs, :nodematcher)
        if !kwargs[:nodematcher](g, h)
            # println("Infeasible: node attribute mismatch")
            return false
        end
    end
    if haskey(kwargs, :edgematcher)
        for (inc, adj) in neighbors(state.G, g)
            haskey(state.g_core, adj) || continue
            hinc = findedgekey(state.H, h, state.g_core[adj])
            if !kwargs[:edgematcher](inc, hinc)
                # println("Infeasible: edge attribute mismatch")
                return false
            end
        end
    end
    return true
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
        v == depth && delete!(state.g_term, k)
    end
    for (k, v) in state.h_term
        v == depth && delete!(state.h_term, k)
    end
end


function yieldallnodes!(state::VF2State)
    push!(state.mappings, copy(state.g_core))
    return
end

function yieldfirstnode!(state::VF2State)
    push!(state.mappings, copy(state.g_core))
    state.status = :done
    return
end

function yieldedgefunc(G, H, yieldtype)
    return function (state)
        delta_y_mismatch(G, H, state.g_core) && return
        push!(state.mappings, copy(state.g_core))
        yieldtype == :first && (state.status = :done)
        return
    end
end


function expand!(state::VF2State; kwargs...)
    state.status == :done && return
    # Recursive
    # println("depth $(length(state.g_core))")
    if length(state.g_core) == nodecount(state.H)
        # println("done $(state.g_core)")
        kwargs[:reportfunc](state)
        return
    elseif state.expire !== nothing && time_ns() > state.expire
        state.status = :timedout
        # println("Timed out")
        return
    end
    for (g, h) in candidatepairs(state; kwargs...)
        # println("candidates $(g) $(h)")
        is_feasible(state, g, h; kwargs...) || continue
        is_semantic_feasible(state, g, h; kwargs...) || continue
        updatestate!(state, g, h)
        # println("g_core $(state.g_core)")
        mp = expand!(state; kwargs...)
        # println("restored $(state.g_core)")
        restore!(state, g, h)
    end
    return
end



function isomorphismitervf2(G::T1, H::T2; timeout=nothing, yield=:all, kwargs...
        ) where {T1<:UndirectedGraph,T2<:UndirectedGraph}
    (nodecount(G) == 0 || nodecount(H) == 0) && return Dict{Int,Int}()
    state = VF2State{T1,T2}(G, H, timeout=timeout)
    yieldfunc = yield == :all ? yieldallnodes! : yieldfirstnode!
    expand!(state, reportfunc=yieldfunc; kwargs...)
    if state.status == :timedout
        throw(ErrorException("Timeout reached"))
    end

    return state.mappings
end


function edgeisomorphismitervf2(
            G::T1, H::T2; timeout=nothing, yield=:all,
            nodematcher=(g,h)->true , edgematcher=(g,h)->true,
            kwargs...
        ) where {T1<:UndirectedGraph,T2<:UndirectedGraph}
    (edgecount(G) == 0 || edgecount(H) == 0) && return Dict{Int,Int}()
    lg = linegraph(G)
    lh = linegraph(H)
    state = VF2State{LineGraph,LineGraph}(lg, lh, timeout=timeout)
    lgematch = lgedgematcher(lg, lh, nodematcher)
    lgnmatch = lgnodematcher(lg, lh, nodematcher, edgematcher)
    yieldfunc = yieldedgefunc(G, H, yield)
    expand!(state, nodematcher=lgnmatch, edgematcher=lgematch,
            reportfunc=yieldfunc; kwargs...)
    if state.status == :timedout
        throw(ErrorException("Timeout reached"))
    end
    return state.mappings
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
    res = isomorphismitervf2(G, H, mode=:Isomorphism, yield=:first; kwargs...)
    return isempty(res) ? nothing : res[1]
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
    res = isomorphismitervf2(G, H, mode=:Subgraph, yield=:first; kwargs...)
    return isempty(res) ? nothing : res[1]
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
The returned iterator has `ig => ih` pairs that correspond to the indices of matching
edges in `G` and `H`, respectively.

See [`MolecularGraph.edgesubgraph`](@ref) to construct the subgraphs that result from the match.
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
    res = edgeisomorphismitervf2(G, H, mode=:Subgraph, yield=:first; kwargs...)
    return isempty(res) ? nothing : res[1]
end

"""
    isedgesubgraphmatch(G::AbstractGraph, H::AbstractGraph; kwargs...) -> Bool

Return true if a node induced subgraph of `G` and `H` are isomorphic.
"""
isedgesubgraphmatch(G, H; kwargs...
    ) = edgesubgraphmatch(G, H; kwargs...) !== nothing
