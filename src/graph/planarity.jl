#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    is_planar,
    is_outerplanar


struct DFSState
    graph::UDGraph # TODO: use type parameter

    rank::Dict{Int,Int} # tree node n, dfsindex(n)
    isthick::Dict{Int,Bool} # edge e, isthick(e)
    inedge::Dict{Int,Int} # tree node n, inedge(n)

    treeedge::Vector{Int} # tree edges in dfs order
    cotree::Set{Int} # cotree edges
    loworder::Dict{Int,Vector{Int}} # tree edge e, list of outedge o in low(o) order
    source::Dict{Int,Int} # tree edge, source node dfsindex
    low::Dict{Int,Int} # cotree edge e, low(e)

    function DFSState(graph)
        new(graph, Dict(), Dict(), Dict(), [], Set(), Dict(), Dict(), Dict())
    end
end


struct DataStructure
    alike::Vector{Int}
    opposite::Vector{Int}
end

DataStructure() = DataStructure([], [])
DataStructure(cotree::Int) = DataStructure([cotree], [])


dfs!(state::DFSState) = dfs!(state, pop!(nodeset(state.graph)))

function dfs!(state::DFSState, u::Int)
    state.rank[u] = length(state.rank) + 1
    buckets = Dict{Int,Vector{Int}}() # low(e), edges
    lows = Dict{Int,Int}() # edge, lownode(e)
    for (v, e) in neighbors(state.graph, u)
        e == get(state.inedge, u, nothing) && continue # Predecessor
        e in state.cotree && continue # Visited
        if !(v in keys(state.rank))
            # Tree node
            push!(state.treeedge, e)
            state.source[e] = state.rank[u]
            state.inedge[v] = e
            low = dfs!(state, v)
        else
            # Cotree node
            low = v
            push!(state.cotree, e)
            state.low[e] = state.rank[v]
            state.isthick[e] = false
        end
        lows[e] = low
        rank = state.rank[low]
        rank in keys(buckets) || (buckets[rank] = Int[])
        push!(buckets[rank], e)
    end
    u in keys(state.inedge) || return # Root
    # low(e) ordering
    inedge = state.inedge[u]
    fringes = [i for i in values(state.low) if i < state.rank[u]]
    state.isthick[inedge] = length(unique(fringes)) > 1
    sorted = Int[]
    for k in sort(collect(keys(buckets)))
        append!(sorted, sort(buckets[k], by=x->state.isthick[x]))
    end
    state.loworder[inedge] = sorted
    return lows[state.loworder[inedge][1]]
end


function monochromatic(a, b, state)
    return isempty(a) || isempty(b) || state.low[a[end]] <= state.low[b[1]]
end


function merge!(ds1::DataStructure, ds2::DataStructure, state::DFSState)
    if isempty(ds1.alike) && isempty(ds1.opposite)
        # Initialize
        append!(ds1.alike, ds2.alike)
        append!(ds1.opposite, ds2.opposite)
        return true
    end
    if monochromatic(ds1.alike, ds2.alike, state)
        append!(ds1.alike, ds2.alike)
        isempty(ds2.opposite) && return true
    elseif monochromatic(ds1.opposite, ds2.alike, state)
        append!(ds1.opposite, ds2.alike)
        isempty(ds2.opposite) && return true
    end
    return false
end


function remove!(ds::DataStructure, state::DFSState, edges)
    setdiff!(ds.alike, edges)
    setdiff!(ds.opposite, edges)
    if (monochromatic(ds.alike, ds.opposite, state)
            || monochromatic(ds.opposite, ds.alike, state))
        append!(ds.alike, ds.opposite)
        empty!(ds.opposite)
    end
end


"""
    is_planar(graph::UDGraph) -> Bool

Return whether the graph is planar.

# Reference

1. de Fraysseix, H., & Ossona de Mendez, P. (2012). Trémaux trees and planarity.
   European Journal of Combinatorics, 33(3), 279–293.
   https://doi.org/10.1016/j.ejc.2011.09.012
"""
function is_planar(graph::UDGraph)
    # Fast filter
    edgecount(graph) < 9 && return true
    edgecount(graph) > nodecount(graph) * 3 - 6 && return false

    # TODO: if not biconnected, decompose it into biconnected components
    g = newgraph(graph) # TODO: does this improve performance?
    # Do DFS to determine treeedge, cotree, loworder, source, low
    state = DFSState(g)
    dfs!(state)
    # println(state)
    # L-R check
    ds = Dict(c => DataStructure(c) for c in state.cotree)
    for e in reverse(state.treeedge)
        ds[e] = DataStructure()
        for i in state.loworder[e]
            # println("stem: $(i), alike: $(ds[i].alike), oppo: $(ds[i].opposite)")
            merge!(ds[e], ds[i], state) || return false
        end
        # println("head: $(n), arrow: $(e), alike: $(ds[e].alike), oppo: $(ds[e].opposite)")
        finished = filter(x -> state.low[x] == state.source[e], state.cotree)
        # println("finished: $(finished)")
        remove!(ds[e], state, finished)
    end
    return true
end


"""
    is_outerplanar(graph::UDGraph) -> Bool

Return whether the graph is outerplanar. The outerplanarity test is based on
a planarity test (see [`is_planar`](@ref)).
"""
function is_outerplanar(graph::UDGraph)
    # Fast filter
    edgecount(graph) < 6 && return true
    edgecount(graph) > nodecount(graph) * 2 - 3 && return false

    g = newgraph(graph)
    nkeys = nodekeys(g)
    newnode = nodecount(g) + 1
    updatenode!(g, nodetype(g)(), newnode)
    ecnt = edgecount(g)
    for n in nkeys
        ecnt += 1
        updateedge!(g, edgetype(g)(newnode, n), ecnt)
    end
    return is_planar(g)
end
