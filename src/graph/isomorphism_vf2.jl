#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

"""
    VF2Matcher{T,U,G<:SimpleGraph{T},H<:SimpleGraph{U}}

Lazy iterator that generate all isomorphism mappings between `G` and `H`.

`matchtype` should be one of the followings
  - `:isomorphic`: G is isomorphic to H
  - `:subgraph_isomorphic`: a node-induced subgraph of G is isomorphic to H
  - `:monomorphic`: a subgraph of G is monomorphic to H

# Options

- `vmatch::Function`: a function for semantic node attribute matching (default: (a, b) -> true)
- `ematch::Function`: a function for semantic edge attribute matching (default: (a, b) -> true)
- `mandatory::Dict{Int,Int}`: mandatory node mapping (if matchtype=:edgeinduced, edge mapping)
- `forbidden::Dict{Int,Int}`: forbidden node mapping (if matchtype=:edgeinduced, edge mapping)
- `timeout::Union{Int,Nothing}`: if specified, abort vf2 calculation when the time reached and return empty iterator (default: 10 seconds)

"""
mutable struct VF2Matcher{T,U,G<:SimpleGraph{T},H<:SimpleGraph{U}}
    g::G
    h::H
    matchtype::Symbol
    vmatch::Function
    ematch::Function
    mandatory::Dict{T,U}
    forbidden::Dict{T,U}
    expire::Union{UInt64,Nothing}  # UInt64, nanoseconds

    stack::Vector{Tuple{Symbol,T,U}}
    currentpair::Union{Tuple{T,U},Nothing}
    g_core::Dict{T,U}
    h_core::Dict{T,U}
    g_term::Dict{T,U}
    h_term::Dict{T,U}
    timedout::Bool

    function VF2Matcher{T,U,G,H}(
                g::G, h::H, matchtype::Symbol;
                vmatch=(a, b)->true, ematch=(a, b)->true,
                mandatory=Dict(), forbidden=Dict(), timeout=nothing
            ) where {T,U,G<:SimpleGraph{T},H<:SimpleGraph{U}}
        expire = (timeout === nothing ?
            nothing : (time_ns() + timeout * 1_000_000_000)::UInt64)
        return new(
            g, h, matchtype, vmatch, ematch,
            mandatory, forbidden, expire,
            [], nothing, Dict(), Dict(), Dict(), Dict(), false
        )
    end
end

VF2Matcher(g::SimpleGraph, h::SimpleGraph, matchtype::Symbol; kwargs...
    ) = VF2Matcher{eltype(g),eltype(h),typeof(g),typeof(h)}(g, h, matchtype; kwargs...)


function candidatepairs(iter::VF2Matcher{T,U,G,H}) where {T,U,G,H}
    # Mandatory pair
    md = setdiff(keys(iter.mandatory), keys(iter.g_core))
    if !isempty(md)
        n = pop!(md)
        return [(n, iter.mandatory[n])]
    end

    pairs = Tuple{T,U}[]
    g_cand = setdiff(keys(iter.g_term), keys(iter.g_core))
    h_cand = setdiff(keys(iter.h_term), keys(iter.h_core))
    if isempty(g_cand) || isempty(h_cand)
        # New connected component
        g_cand = setdiff(Set(vertices(iter.g)), keys(iter.g_core))
        h_cand = setdiff(Set(vertices(iter.h)), keys(iter.h_core))
    end
    if !isempty(h_cand)
        # improved pivot for monomorphism match with many isolated nodes
        h_sorted = sort!(collect(h_cand))
        h_min = h_sorted[argmax([degree(iter.h, n) for n in h_sorted])]
        for i in g_cand
            # Forbidden pair
            if haskey(iter.forbidden, i) && iter.forbidden[i] == h_min
                continue
            end
            push!(pairs, (i, h_min))
        end
    end
    return pairs
end


function is_feasible(iter::VF2Matcher, gv, hv)
    # Note: assume simple graph (no self loops and multi edges)
    g_nbrs = neighbors(iter.g, gv)
    h_nbrs = neighbors(iter.h, hv)
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

        g_new_count = length(setdiff(Set(vertices(iter.g)), keys(iter.g_term)))
        h_new_count = length(setdiff(Set(vertices(iter.h)), keys(iter.h_term)))
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


function is_semantic_feasible(iter::VF2Matcher, gv, hv)
    if !iter.vmatch(gv, hv)
        @debug "Infeasible: node attribute mismatch"
        return false
    end
    for nbr in neighbors(iter.g, gv)
        haskey(iter.g_core, nbr) || continue
        (iter.matchtype === :monomorphic
            && !has_edge(iter.h, hv, iter.g_core[nbr])) && continue
        if !iter.ematch(
                u_edge(iter.g, gv, nbr),
                u_edge(iter.h, hv, iter.g_core[nbr]))
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
    g_nbrset = union([neighbors(iter.g, n) for n in keys(iter.g_core)]...)
    for n in setdiff(g_nbrset, keys(iter.g_term))
        iter.g_term[n] = depth
    end
    h_nbrset = union([neighbors(iter.h, n) for n in keys(iter.h_core)]...)
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
        push!(iter.stack, (:expand, g, h))
    end
    while !isempty(iter.stack)
        command, g1, h1 = pop!(iter.stack)
        if command === :expand
            is_feasible(iter, g1, h1) || continue
            is_semantic_feasible(iter, g1, h1) || continue
        elseif command === :restore
            @debug "restored" iter.g_core
            restore!(iter, g1, h1)
            continue
        end
        @debug "g_core" iter.g_core
        expand!(iter, g1, h1)
        iter.currentpair = (g1, h1)
        @debug "depth" length(iter.g_core)
        if length(iter.g_core) == nv(iter.h)
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
            push!(iter.stack, (:expand, g, h))
        end
    end
    return  # no match found
end

Base.IteratorSize(::Type{<:VF2Matcher}) = Base.SizeUnknown()
Base.IteratorEltype(::Type{<:VF2Matcher}) = Base.EltypeUnknown()


# Graph isomorphism

"""
    isomorphisms(g::SimpleGraph, h::SimpleGraph; kwargs...) -> Iterator

Return an iterator that generate isomorphic mappings between `G` and `H`.
The returned iterator has `ig => ih` pairs that correspond to the indices of matching nodes in `G` and `H`, respectively.
"""
isomorphisms(g, h; kwargs...) = VF2Matcher(g, h, :isomorphic; kwargs...) 


"""
    is_isomorphic(g::SimpleGraph, h::SimpleGraph; kwargs...) -> Bool

Return whether `G` and `H` are isomorphic.
"""
is_isomorphic(g, h; kwargs...) = !isempty(isomorphisms(g, h; kwargs...))


# Node-induced subgraph isomorphism

"""
    nodesubgraph_isomorphisms(g::SimpleGraph, h::SimpleGraph; kwargs...) -> Iterator

Return an iterator that generate isomorphic mappings between `H` and node-induced subgraphs of `G`.
The returned iterator has `ig => ih` pairs that correspond to the indices of matching nodes in `G` and `H`, respectively.
"""
nodesubgraph_isomorphisms(g, h; kwargs...
    ) = VF2Matcher(g, h, :subgraph_isomorphic; kwargs...) 


"""
    nodesubgraph_is_isomorphic(g::SimpleGraph, h::SimpleGraph; kwargs...) -> Bool

Return whether a node-induced subgraph of `G` is isomorphic to `H`.
"""
nodesubgraph_is_isomorphic(g, h; kwargs...
    ) = !isempty(nodesubgraph_isomorphisms(g, h; kwargs...))



# Edge-induced subgraph isomorphism

"""
    edgesubgraph_isomorphisms(
        g::SimpleGraph, h::SimpleGraph;
        vmatch=(g,h)->true, ematch=(g,h)->true,
        kwargs...) -> Iterator

Return an iterator that generate isomorphic mappings between `H` and edge-induced subgraphs of `G`.
The returned iterator has `ig => ih` pairs that correspond to the indices of matching edges
in `G` and `H`, respectively.

`vmatch` and `ematch` control the features needed to be counted as a match.

See `Graphs.induced_subgraph` to construct the subgraphs that result from the match.
"""
function edgesubgraph_isomorphisms(
        g, h; vmatch=(gv,hv)->true, ematch=(ge,he)->true,
        ggmatch=(n1,n2)->true, hhmatch=(n1,n2)->true, kwargs...)
    lg, grev, gsh = line_graph(g)
    lh, hrev, hsh = line_graph(h)
    lgnode = lgvmatch(lg, lh, grev, hrev, vmatch, ematch)
    lgedge = lgematch(lg, lh, gsh, hsh, vmatch)
    gd, gy, hd, hy = (
        delta_edges(g, ggmatch), y_edges(g, ggmatch), delta_edges(h, hhmatch), y_edges(h, hhmatch))
    matcher = VF2Matcher(
        lg, lh, :subgraph_isomorphic;
        vmatch=lgnode, ematch=lgedge, kwargs...)
    revmaped = Iterators.map(matcher) do mapping
        return Dict(grev[m] => hrev[n] for (m, n) in mapping)
    end
    return Iterators.filter(revmaped) do mapping
        return delta_y_test(mapping, gd, gy, hd, hy)
    end
end


"""
    edgesubgraph_is_isomorphic(g::SimpleGraph, h::SimpleGraph; kwargs...) -> Bool

Return whether an edge-induced subgraph of `G` is isomorphic to `H`.
"""
edgesubgraph_is_isomorphic(g, h; kwargs...) = !isempty(edgesubgraph_isomorphisms(g, h; kwargs...))


# Subgraph monomorphism

"""
    subgraph_monomorphisms(g::SimpleGraph, h::SimpleGraph; kwargs...) -> Iterator

Generate monomorphism mappings between `H` and subgraphs of `G`.
The returned iterator has `ig => ih` pairs that correspond to the indices of
matching nodes in `G` and `H`, respectively.
"""
subgraph_monomorphisms(g, h; kwargs...) = VF2Matcher(g, h, :monomorphic; kwargs...)


"""
    subgraph_is_monomorphic(g::SimpleGraph, h::SimpleGraph; kwargs...) -> Bool

Return whether a subgraph of `G` is monomorphic to `H`.
"""
subgraph_is_monomorphic(g, h; kwargs...) = !isempty(subgraph_monomorphisms(g, h; kwargs...))
