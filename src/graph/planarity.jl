#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    planaritytest, outerplanaritytest,
    isplanar, isouterplanar


PlanarityTestDS{T} = Vector{Vector{Vector{T}}}


struct PlanarityTestState{T}
    graph::SimpleGraph{T}

    rank::Dict{T,Int} # tree node n, dfsindex(n)
    inedge::Dict{T,Edge{T}} # tree node n, inedge(n)
    isthick::Dict{Edge{T},Bool} # edge e, isthick(e)

    source::Dict{Edge{T},Int} # edge e, source node dfsindex
    loworder::Dict{Edge{T},Vector{Edge{T}}} # tree edge e, list of outedge o in low(o) order
    cotree::Dict{Edge{T},Int} # cotree edge e, low(e)
    treeedge::Vector{Edge{T}} # tree edges in dfs order
end

PlanarityTestState(g::SimpleGraph{T}
    ) where T = PlanarityTestState{T}(g, Dict(), Dict(), Dict(), Dict(), Dict(), Dict(), [])


function dfs!(state::PlanarityTestState{T}, u::Int) where T
    state.rank[u] = length(state.rank) + 1
    buckets = Dict{T,Vector{Edge{T}}}() # low(e), edges
    lows = Dict{Edge{T},T}() # edge, lownode(e)
    for nbr in neighbors(state.graph, u)
        inc = u_edge(state.graph, u, nbr)
        inc == get(state.inedge, u, nothing) && continue # Predecessor
        haskey(state.cotree, inc) && continue # Visited
        state.source[inc] = state.rank[u]
        if !haskey(state.rank, nbr)
            # Tree node
            push!(state.treeedge, inc)
            state.inedge[nbr] = inc
            low = dfs!(state, nbr)
        else
            # Cotree node
            low = nbr
            state.cotree[inc] = state.rank[nbr]
            state.isthick[inc] = false
        end
        low === nothing && continue
        lows[inc] = low
        rank = state.rank[low]
        haskey(buckets, rank) || (buckets[rank] = Edge{T}[])
        push!(buckets[rank], inc)
    end
    u in keys(state.inedge) || return # Root
    # low(e) ordering
    inedge = state.inedge[u]
    fringes = [i for i in values(state.cotree) if i < state.rank[u]]
    state.isthick[inedge] = length(unique(fringes)) > 1
    sorted = Edge{T}[]
    for k in sort(collect(keys(buckets)))
        append!(sorted, sort(buckets[k], by=x->state.isthick[x]))
    end
    state.loworder[inedge] = sorted
    isempty(lows) && return # No cycles
    return lows[state.loworder[inedge][1]]
end


function merge_ds!(ds1::PlanarityTestDS{T}, ds2::PlanarityTestDS{T}, cotree) where T
    if isempty(ds1)
        append!(ds1, ds2)
        @debug "Trunk"
        return true
    end
    @assert !isempty(ds2)
    b1 = cotree[ds1[1][1][1]]
    b2 = cotree[ds2[1][1][1]]
    newds = PlanarityTestDS{T}()
    lowtop = 0
    while lowtop <= b2 && !isempty(ds1)
        pair = popfirst!(ds1)
        topl = isempty(pair[1]) ? 0 : cotree[pair[1][end]]
        topr = isempty(pair[2]) ? 0 : cotree[pair[2][end]]
        if topl > b2 && topr > b2
            @debug "Not planer: j+1 case"
            return false
        end
        lowtop = max(topl, topr)
        @assert lowtop != 0
        newl = topl >= topr ? 1 : 2
        newr = topl >= topr ? 2 : 1
        push!(newds, [pair[newl], pair[newr]])
    end
    # newds[end] is now j + 1
    if lowtop <= b2
        push!(newds, [[], []])
    end

    fused1 = T[]
    while !isempty(ds1)
        pair = popfirst!(ds1)
        if !isempty(pair[1]) && !isempty(pair[2])
            @debug "Not planer: should be trichromatic"
            return false
        end
        append!(fused1, pair[1], pair[2])
    end
    append!(newds[end][1], fused1)

    fused2 = T[]
    for (i, pair) in enumerate(ds2)
        if !isempty(pair[1]) && !isempty(pair[2])
            @debug "Not planer: dichromatic branch cannot be merged"
            return false
        end
        cell = i == 1 && b1 == b2 ? newds[1][1] : fused2
        append!(cell, pair[1], pair[2])
    end
    whichtofuse = isempty(newds[end][1]) ? 1 : 2
    append!(newds[end][whichtofuse], fused2)
    if isempty(newds[end][1]) && isempty(newds[end][2])
        # ds2 has only a b1 == b2 edge
        pop!(newds)
    end
    # Already ds1 has used up
    append!(ds1, newds)
    @debug "Merged: ds $(ds1)"
    return true
end


function planaritytest(g::SimpleGraph{T}) where T
    # Do DFS to determine treeedge, cotree, loworder, source
    state = PlanarityTestState(g)
    nodes = Set(vertices(g))
    while !isempty(nodes)
        n = pop!(nodes)
        dfs!(state, n)
        setdiff!(nodes, keys(state.rank))
    end
    @debug state

    # L-R check
    ds = Dict(c => Vector{Vector{Edge{T}}}[[[c], []]] for c in keys(state.cotree))
    while !isempty(state.treeedge)
        e = pop!(state.treeedge)
        ds[e] = PlanarityTestDS{Edge{T}}()
        @debug "inedge: $(e)"
        for i in state.loworder[e]
            @debug "stem: $(i), ds: $(ds[i])"
            if !merge_ds!(ds[e], ds[i], state.cotree)
                return false
            end
        end
        finished = filter(x -> state.cotree[x] == state.source[e], keys(state.cotree))
        @debug "removed: $(finished)"
        if !isempty(ds[e])
            setdiff!(ds[e][end][1], finished)
            setdiff!(ds[e][end][2], finished)
            if isempty(ds[e][end][1]) && isempty(ds[e][end][2])
                pop!(ds[e])
            end
        end
    end
    return true
end


function outerplanaritytest(g::SimpleGraph)
    # Add a node connects all other nodes.
    # If the modified graph is planer, the original graph is outerplaner.
    newg = copy(g)
    add_vertex!(newg)
    newnode = vertices(newg)[end]
    for i in vertices(g)
        add_edge!(newg, newnode, i)
    end
    return planaritytest(newg)
end


"""
    isplanar(graph::SimpleGraph) -> Bool

Return whether the graph is planar.

# Reference

1. de Fraysseix, H., & Ossona de Mendez, P. (2012). Trémaux trees and planarity.
   European Journal of Combinatorics, 33(3), 279–293.
   https://doi.org/10.1016/j.ejc.2011.09.012
1. Trémaux trees and planarity. https://arxiv.org/abs/math/0610935
"""
function isplanar(g::SimpleGraph)
    ne(g) < 9 && return true
    ne(g) > nv(g) * 3 - 6 && return false
    return planaritytest(g)
end


"""
    isouterplanar(graph::SimpleGraph) -> Bool

Return whether the graph is outerplanar. The outerplanarity test is based on
a planarity test (see [`isplanar`](@ref)).
"""
function isouterplanar(g::SimpleGraph)
    ne(g) < 6 && return true
    ne(g) > nv(g) * 2 - 6 && return false
    return outerplanaritytest(g)
end
