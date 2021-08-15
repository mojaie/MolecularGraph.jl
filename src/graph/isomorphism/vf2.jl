#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    VF2Matcher,
    isomorphisms, is_isomorphic,
    nodesubgraph_isomorphisms, nodesubgraph_is_isomorphic,
    edgesubgraph_isomorphisms, edgesubgraph_is_isomorphic,
    subgraph_monomorphisms, subgraph_is_monomorphic


"""
    VF2Matcher{T1<:UndirectedGraph,T2<:UndirectedGraph}(
        G::T1, H::T2, matchtype::Symbol; kwargs...
    )

Lazy iterator that generate all isomorphism mappings between `G` and `H`.

`matchtype` should be one of the followings
  - `:isomorphic`: G is isomorphic to H
  - `:subgraph_isomorphic`: a node-induced subgraph of G is isomorphic to H
  - `:monomorphic`: a subgraph of G is monomorphic to H

# Options

- `nodematcher::Function`: a function for semantic node attribute matching (default: (a, b) -> true)
- `edgematcher::Function`: a function for semantic edge attribute matching (default: (a, b) -> true)
- `mandatory::Dict{Int,Int}`: mandatory node mapping (if matchtype=:edgeinduced, edge mapping)
- `forbidden::Dict{Int,Int}`: forbidden node mapping (if matchtype=:edgeinduced, edge mapping)
- `timeout::Union{Int,Nothing}`: if specified, abort vf2 calculation when the time reached and return empty iterator (default: 10 seconds)

"""
mutable struct VF2Matcher{T1<:UndirectedGraph,T2<:UndirectedGraph}
    G::T1
    H::T2
    matchtype::Symbol
    nodematcher::Function
    edgematcher::Function
    mandatory::Dict{Int,Int}
    forbidden::Dict{Int,Int}
    expire::Union{UInt64,Nothing} # UInt64, nanoseconds

    stack::Vector{Tuple{Symbol,Int,Int}}
    currentpair::Union{Tuple{Int,Int},Nothing}
    g_core::Dict{Int,Int}
    h_core::Dict{Int,Int}
    g_term::Dict{Int,Int}
    h_term::Dict{Int,Int}
    timedout::Bool

    function VF2Matcher{T1,T2}(
                G::T1, H::T2, matchtype::Symbol;
                nodematcher=(a, b)->true, edgematcher=(a, b)->true,
                mandatory=Dict(), forbidden=Dict(), timeout=nothing
            ) where {T1<:UndirectedGraph,T2<:UndirectedGraph}
        expire = (timeout === nothing ?
            nothing : (time_ns() + timeout * 1_000_000_000)::UInt64)
        return new(
            G, H, matchtype, nodematcher, edgematcher,
            mandatory, forbidden, expire,
            [], nothing, Dict(), Dict(), Dict(), Dict(), false
        )
    end
end
VF2Matcher(
    G::UndirectedGraph, H::UndirectedGraph, matchtype::Symbol; kwargs...
) = VF2Matcher{typeof(G),typeof(H)}(G, H, matchtype; kwargs...)


function candidatepairs(iter::VF2Matcher)
    # Mandatory pair
    md = setdiff(keys(iter.mandatory), keys(iter.g_core))
    if !isempty(md)
        n = pop!(md)
        return [(n, iter.mandatory[n])]
    end

    pairs = Tuple{Int,Int}[]
    g_cand = setdiff(keys(iter.g_term), keys(iter.g_core))
    h_cand = setdiff(keys(iter.h_term), keys(iter.h_core))
    if isempty(g_cand) || isempty(h_cand)
        # New connected component
        g_cand = setdiff(nodeset(iter.G), keys(iter.g_core))
        h_cand = setdiff(nodeset(iter.H), keys(iter.h_core))
    end
    if !isempty(h_cand)
        # improved pivot for monomorphism match with many isolated nodes
        h_sorted = sort!(collect(h_cand))
        h_min = h_sorted[argmax([degree(iter.H, n) for n in h_sorted])]
        for g in g_cand
            # Forbidden pair
            if haskey(iter.forbidden, g) && iter.forbidden[g] == h_min
                continue
            end
            push!(pairs, (g, h_min))
        end
    end
    return pairs
end


function is_feasible(iter::VF2Matcher, g, h)
    # Note: assume simple graph (no self loops and multi edges)
    g_nbrs = adjacencies(iter.G, g)
    h_nbrs = adjacencies(iter.H, h)
    if iter.matchtype !== :monomorphic
        for n in intersect(g_nbrs, keys(iter.g_core))
            if !(iter.g_core[n] in h_nbrs)
                @debug "Infeasible: neighbor connectivity"
                return false
            end
        end
    end
    for n in intersect(h_nbrs, keys(iter.h_core))
        if !(iter.h_core[n] in g_nbrs)
            @debug "Infeasible: neighbor connectivity"
            return false
        end
    end
    
    if iter.matchtype !== :monomorphic
        g_term_count = length(setdiff(keys(iter.g_term), keys(iter.g_core)))
        h_term_count = length(setdiff(keys(iter.h_term), keys(iter.h_core)))
        if iter.matchtype === :exact && g_term_count != h_term_count
            @debug "Infeasible: terminal set size"
            return false
        elseif iter.matchtype === :subgraph_isomorphic && g_term_count < h_term_count
            @debug "Infeasible: terminal set size"
            return false
        end

        g_new_count = length(setdiff(nodeset(iter.G), keys(iter.g_term)))
        h_new_count = length(setdiff(nodeset(iter.H), keys(iter.h_term)))
        if iter.matchtype === :isomorphic && g_new_count != h_new_count
            @debug "Infeasible: yet unexplored nodes"
            return false
        elseif iter.matchtype === :subgraph_isomorphic && g_new_count < h_new_count
            @debug "Infeasible: yet unexplored nodes"
            return false
        end
    end
    @debug "Syntactically feasible"
    return true
end


function is_semantic_feasible(iter::VF2Matcher, g, h)
    if !iter.nodematcher(g, h)
        @debug "Infeasible: node attribute mismatch"
        return false
    end
    for (inc, adj) in neighbors(iter.G, g)
        haskey(iter.g_core, adj) || continue
        hinc = findedgekey(iter.H, h, iter.g_core[adj])
        iter.matchtype === :monomorphic && hinc === nothing && continue
        if !iter.edgematcher(inc, hinc)
            @debug "Infeasible: edge attribute mismatch"
            return false
        end
    end
    return true
end


function expand!(iter::VF2Matcher, g, h)
    iter.g_core[g] = h
    iter.h_core[h] = g
    depth = length(iter.g_core)
    if !haskey(iter.g_term, g)
        iter.g_term[g] = depth
    end
    if !haskey(iter.h_term, h)
        iter.h_term[h] = depth
    end
    g_nbrset = union([adjacencies(iter.G, n) for n in keys(iter.g_core)]...)
    for n in setdiff(g_nbrset, keys(iter.g_term))
        iter.g_term[n] = depth
    end
    h_nbrset = union([adjacencies(iter.H, n) for n in keys(iter.h_core)]...)
    for n in setdiff(h_nbrset, keys(iter.h_term))
        iter.h_term[n] = depth
    end
    return
end


function restore!(iter::VF2Matcher, g, h)
    depth = length(iter.g_core)
    if g !== nothing && h !== nothing
        delete!(iter.g_core, g)
        delete!(iter.h_core, h)
    end
    for (k, v) in iter.g_term
        v == depth && delete!(iter.g_term, k)
    end
    for (k, v) in iter.h_term
        v == depth && delete!(iter.h_term, k)
    end
end


function Base.iterate(iter::VF2Matcher, state=nothing)
    iter.timedout && return
    if iter.currentpair !== nothing
        g1, h1 = iter.currentpair
        push!(iter.stack, (:restore, g1, h1))
    end
    for (g, h) in candidatepairs(iter)
        @debug "candidates" g h
        is_feasible(iter, g, h) || continue
        is_semantic_feasible(iter, g, h) || continue
        push!(iter.stack, (:expand, g, h))
    end
    while !isempty(iter.stack)
        command, g1, h1 = pop!(iter.stack)
        if command === :restore
            @debug "restored" iter.g_core
            restore!(iter, g1, h1)
            continue
        end
        @debug "g_core" iter.g_core
        expand!(iter, g1, h1)
        iter.currentpair = (g1, h1)
        @debug "depth" length(iter.g_core)
        if length(iter.g_core) == nodecount(iter.H)
            @debug "done" iter.g_core
            return (copy(iter.g_core), state)
        elseif iter.expire !== nothing && time_ns() > iter.expire
            iter.timedout = true
            @debug "Timed out"
            return
        end
        push!(iter.stack, (:restore, g1, h1))
        for (g, h) in candidatepairs(iter)
            @debug "candidates" g h
            is_feasible(iter, g, h) || continue
            is_semantic_feasible(iter, g, h) || continue
            push!(iter.stack, (:expand, g, h))
        end
    end
    return  # no match found
end

# TODO: get rid of generics
Base.IteratorSize(::Type{VF2Matcher{T1,T2}}) where {T1<:UndirectedGraph,T2<:UndirectedGraph} = Base.SizeUnknown()
Base.IteratorEltype(::Type{VF2Matcher{T1,T2}}) where {T1<:UndirectedGraph,T2<:UndirectedGraph} = Base.EltypeUnknown()



# Graph isomorphism

"""
    isomorphisms(G::UndirectedGraph, H::UndirectedGraph; kwargs...) -> Iterator

Return an iterator that generate isomorphic mappings between `G` and `H`.
The returned iterator has `ig => ih` pairs that correspond to the indices of matching nodes in `G` and `H`, respectively.
"""
isomorphisms(G, H; kwargs...) = VF2Matcher(G, H, :isomorphic; kwargs...) 


"""
    is_isomorphic(G::UndirectedGraph, H::UndirectedGraph; kwargs...) -> Bool

Return whether `G` and `H` are isomorphic.
"""
is_isomorphic(G, H; kwargs...) = !isempty(isomorphisms(G, H; kwargs...))


# Node-induced subgraph isomorphism

"""
    nodesubgraph_isomorphisms(G::UndirectedGraph, H::UndirectedGraph; kwargs...) -> Iterator

Return an iterator that generate isomorphic mappings between `H` and node-induced subgraphs of `G`.
The returned iterator has `ig => ih` pairs that correspond to the indices of matching nodes in `G` and `H`, respectively.
"""
nodesubgraph_isomorphisms(G, H; kwargs...
    ) = VF2Matcher(G, H, :subgraph_isomorphic; kwargs...) 


"""
    nodesubgraph_is_isomorphic(G::UndirectedGraph, H::UndirectedGraph; kwargs...) -> Bool

Return whether a node-induced subgraph of `G` is isomorphic to `H`.
"""
nodesubgraph_is_isomorphic(G, H; kwargs...
    ) = !isempty(nodesubgraph_isomorphisms(G, H; kwargs...))



# Edge-induced subgraph isomorphism

"""
    edgesubgraph_isomorphisms(
        G::UndirectedGraph, H::UndirectedGraph;
        nodematcher=(g,h)->true, edgematcher=(g,h)->true,
        kwargs...) -> Iterator

Return an iterator that generate isomorphic mappings between `H` and edge-induced subgraphs of `G`.
The returned iterator has `ig => ih` pairs that correspond to the indices of matching edges in `G` and `H`, respectively.

`nodematcher` and `edgematcher` control the features needed to be counted as a match.

See [`Graph.edgesubgraph`](@ref) to construct the subgraphs that result from the match.
"""
function edgesubgraph_isomorphisms(
        G, H; nodematcher=(g,h)->true, edgematcher=(g,h)->true, kwargs...)
    lg = linegraph(G)
    lh = linegraph(H)
    lgedge = lgedgematcher(lg, lh, nodematcher)
    lgnode = lgnodematcher(lg, lh, nodematcher, edgematcher)
    matcher = VF2Matcher(
        lg, lh, :subgraph_isomorphic;
        nodematcher=lgnode, edgematcher=lgedge, kwargs...)
    return Iterators.filter(matcher) do mapping
        return !delta_y_mismatch(G, H, mapping)
    end
end


"""
    edgesubgraph_is_isomorphic(
        G::UndirectedGraph, H::UndirectedGraph; kwargs...) -> Bool

Return whether an edge-induced subgraph of `G` is isomorphic to `H`.
"""
edgesubgraph_is_isomorphic(G, H; kwargs...
    ) = !isempty(edgesubgraph_isomorphisms(G, H; kwargs...))


# Subgraph monomorphism

"""
    subgraph_monomorphisms(
        G::UndirectedGraph, H::UndirectedGraph; kwargs...) -> Iterator

Generate monomorphism mappings between `H` and subgraphs of `G`.
The returned iterator has `ig => ih` pairs that correspond to the indices of matching nodes in `G` and `H`, respectively.
"""
subgraph_monomorphisms(G, H; kwargs...) = VF2Matcher(G, H, :monomorphic; kwargs...)


"""
    subgraph_is_monomorphic(
        G::UndirectedGraph, H::UndirectedGraph; kwargs...) -> Bool

Return whether a subgraph of `G` is monomorphic to `H`.
"""
subgraph_is_monomorphic(G, H; kwargs...
    ) = !isempty(subgraph_monomorphisms(G, H; kwargs...))
