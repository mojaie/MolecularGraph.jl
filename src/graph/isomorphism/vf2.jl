#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    VF2Matcher,
    allexactmatches, exactmatch, isexactmatch,
    allsubgraphmatches, subgraphmatch, issubgraphmatch,
    alledgesubgraphmatches, edgesubgraphmatch, isedgesubgraphmatch


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
                G::T1, H::T2, matchtype;
                nodematcher=(a, b)->true, edgematcher=(a, b)->true,
                mandatory=nothing, forbidden=nothing, timeout=nothing
            ) where {T1<:UndirectedGraph,T2<:UndirectedGraph}
        md = mandatory === nothing ? Dict() : mandatory
        fb = forbidden === nothing ? Dict() : forbidden
        expire = (timeout === nothing ?
            nothing : (time_ns() + timeout * 1_000_000_000)::UInt64)
        return new(
            G, H, matchtype, nodematcher, edgematcher, md, fb, expire,
            [], nothing, Dict(), Dict(), Dict(), Dict(), false
        )
    end
end
VF2Matcher(G::UndirectedGraph, H::UndirectedGraph; kwargs...) = VF2Matcher{typeof(G),typeof(H)}(G, H; kwargs...)


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
        h_min = minimum(h_cand)
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
    # TODO: assume no self loop
    g_nbrs = adjacencies(iter.G, g)
    h_nbrs = adjacencies(iter.H, h)
    for n in intersect(g_nbrs, keys(iter.g_core))
        if !(iter.g_core[n] in h_nbrs)
            @debug "Infeasible: neighbor connectivity"
            return false
        end
    end
    for n in intersect(h_nbrs, keys(iter.h_core))
        if !(iter.h_core[n] in g_nbrs)
            @debug "Infeasible: neighbor connectivity"
            return false
        end
    end
    g_term_count = length(setdiff(keys(iter.g_term), keys(iter.g_core)))
    h_term_count = length(setdiff(keys(iter.h_term), keys(iter.h_core)))
    if iter.matchtype === :exact && g_term_count != h_term_count
        @debug "Infeasible: terminal set size"
        return false
    elseif iter.matchtype === :subgraph && g_term_count < h_term_count
        @debug "Infeasible: terminal set size"
        return false
    end
    g_new_count = length(setdiff(nodeset(iter.G), keys(iter.g_term)))
    h_new_count = length(setdiff(nodeset(iter.H), keys(iter.h_term)))
    if iter.matchtype === :exact && g_new_count != h_new_count
        @debug "Infeasible: yet unexplored nodes"
        return false
    elseif iter.matchtype === :subgraph && g_new_count < h_new_count
        @debug "Infeasible: yet unexplored nodes"
        return false
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

Base.IteratorSize(::Type{VF2Matcher}) = Base.SizeUnknown()
Base.IteratorEltype(::Type{VF2Matcher}) = Base.EltypeUnknown()



# Graph isomorphism

"""
    allexactmatches(G::AbstractGraph, H::AbstractGraph; kwargs...) -> Iterator

Generate isomorphism mappings between `G` and `H`. If no match found, return
nothing.
"""
allexactmatches(G::T1, H::T2; kwargs...
    ) where {T1<:UndirectedGraph,T2<:UndirectedGraph
    } = VF2Matcher{T1,T2}(G, H, :exact; kwargs...) 

"""
    exactmatch(G::AbstractGraph, H::AbstractGraph; kwargs...) -> Dict{Int,Int}

Return an isomorphism mapping between `G` and `H`. If no match found, return
nothing.
"""
exactmatch(G, H; kwargs...) = iterate(allexactmatches(G, H; kwargs...))

"""
    isexactmatch(G::AbstractGraph, H::AbstractGraph; kwargs...) -> Bool

Return true if `G` and `H` are isomorphic.
"""
isexactmatch(G, H; kwargs...) = exactmatch(G, H; kwargs...) !== nothing



"""
    subgraphmatches(G::AbstractGraph, H::AbstractGraph; kwargs...) -> Iterator

Generate node induced subgraph isomorphism mappings between `G` and `H`.

# Keyword arguments

- nodematcher(Function): node matcher function that takes two node indices as
arguments.
- edgematcher(Function): edge matcher function that takes two edge indices as
arguments.
- mandatory(Dict{Int,Int}): mandatory node matches (available for only VF2)
- forbidden(Dict{Int,Int}):
    forbidden node matches (available for only VF2)
"""
allsubgraphmatches(G::T1, H::T2; kwargs...
    ) where {T1<:UndirectedGraph,T2<:UndirectedGraph
    } = VF2Matcher{T1,T2}(G, H, :subgraph; kwargs...) 

"""
    subgraphmatch(
        G::AbstractGraph, H::AbstractGraph; kwargs...) -> Dict{Int,Int}

Return a subgraph isomorphism mapping between `G` and `H`. If no match found,
return nothing.
"""
subgraphmatch(G, H; kwargs...) = iterate(allsubgraphmatches(G, H; kwargs...))

"""
    issubgraphmatch(G::AbstractGraph, H::AbstractGraph; kwargs...) -> Bool

Return true if `H` and a node induced subgraph of `G` are isomorphic.
"""
issubgraphmatch(G, H; kwargs...) = subgraphmatch(G, H; kwargs...) !== nothing



# Edge induced subgraph isomorphism

"""
    alledgesubgraphmatches(
        G::AbstractGraph, H::AbstractGraph;
        nodematcher=(g,h)->true , edgematcher=(g,h)->true,
        kwargs...) -> Iterator

Generate edge induced subgraph isomorphism mappings between `G` and `H`.
The returned iterator has `ig => ih` pairs that correspond to the indices of matching
edges in `G` and `H`, respectively.

`nodematcher` and `edgematcher` control the features needed to be counted as a match.
[`atommatch`](@ref) and [`bondmatch`](@ref) can be used to construct these functions.

See [`MolecularGraph.edgesubgraph`](@ref) to construct the subgraphs that result from the match.
"""
function alledgesubgraphmatches(
            G::T1, H::T2,
            nodematcher=(g,h)->true, edgematcher=(g,h)->true; kwargs...
        ) where {T1<:UndirectedGraph,T2<:UndirectedGraph}
    lg = linegraph(G)
    lh = linegraph(H)
    lgedge = lgedgematcher(lg, lh, nodematcher)
    lgnode = lgnodematcher(lg, lh, nodematcher, edgematcher)
    matcher = VF2Matcher{LineGraph,LineGraph}(
        lg, lh, :subgraph, nodematcher=lgnode, edgematcher=lgedge; kwargs...
    )
    return Iterators.filter(matcher) do mapping
        return !delta_y_mismatch(G, H, mapping)
    end
end

"""
    edgesubgraphmatch(
        G::AbstractGraph, H::AbstractGraph; kwargs...) -> Dict{Int,Int}

Return a edge induced subgraph isomorphism mapping between `G` and `H`.
If no match found, return nothing.
"""
edgesubgraphmatch(G, H; kwargs...
    ) = iterate(alledgesubgraphmatches(G, H; kwargs...))

"""
    isedgesubgraphmatch(G::AbstractGraph, H::AbstractGraph; kwargs...) -> Bool

Return true if a node induced subgraph of `G` and `H` are isomorphic.
"""
isedgesubgraphmatch(G, H; kwargs...
    ) = edgesubgraphmatch(G, H; kwargs...) !== nothing
