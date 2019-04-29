#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    isplanar,
    isouterplanar


DataStructure = Vector{Vector{Vector{Int}}}


struct DFSState
    graph::PlainGraph

    rank::Dict{Int,Int} # tree node n, dfsindex(n)
    inedge::Dict{Int,Int} # tree node n, inedge(n)
    isthick::Dict{Int,Bool} # edge e, isthick(e)

    source::Dict{Int,Int} # edge e, source node dfsindex
    loworder::Dict{Int,Vector{Int}} # tree edge e, list of outedge o in low(o) order
    cotree::Dict{Int,Int} # cotree edge e, low(e)
    treeedge::Vector{Int} # tree edges in dfs order

    function DFSState(graph)
        new(graph, Dict(), Dict(), Dict(), Dict(), Dict(), Dict(), [])
    end
end


dfs!(state::DFSState) = dfs!(state, 1)

function dfs!(state::DFSState, u::Int)
    state.rank[u] = length(state.rank) + 1
    buckets = Dict{Int,Vector{Int}}() # low(e), edges
    lows = Dict{Int,Int}() # edge, lownode(e)
    for (inc, adj) in neighbors(state.graph, u)
        inc == get(state.inedge, u, nothing) && continue # Predecessor
        haskey(state.cotree, inc) && continue # Visited
        state.source[inc] = state.rank[u]
        if !haskey(state.rank, adj)
            # Tree node
            push!(state.treeedge, inc)
            state.inedge[adj] = inc
            low = dfs!(state, adj)
        else
            # Cotree node
            low = adj
            state.cotree[inc] = state.rank[adj]
            state.isthick[inc] = false
        end
        lows[inc] = low
        rank = state.rank[low]
        haskey(buckets, rank) || (buckets[rank] = Int[])
        push!(buckets[rank], inc)
    end
    u == 1 && return # Root
    # low(e) ordering
    inedge = state.inedge[u]
    fringes = [i for i in values(state.cotree) if i < state.rank[u]]
    state.isthick[inedge] = length(unique(fringes)) > 1
    sorted = Int[]
    for k in sort(collect(keys(buckets)))
        append!(sorted, sort(buckets[k], by=x->state.isthick[x]))
    end
    state.loworder[inedge] = sorted
    return lows[state.loworder[inedge][1]]
end


function merge!(ds1::DataStructure, ds2::DataStructure,
        cotree::Dict{Int,Int}; verbose=false)
    if isempty(ds1)
        append!(ds1, ds2)
        verbose && println("Trunk cells")
        return true
    end
    monoc(a, b) = cotree[a[end]] <= cotree[b[1]]
    branch = Int[]
    rootcells = DataStructure()
    # Collect branch cells
    for (alike, oppo) in ds2
        if !isempty(oppo)
            verbose && println("Not planer: dichromatic branch found")
            return false
        end
        if length(alike) == 1 && cotree[ds1[1][1][1]] == cotree[alike[1]]
            c = Vector{Int}[alike, []]
            verbose && println("Root cell $(c)")
            push!(rootcells, c)
        else
            append!(branch, alike)
        end
    end
    if isempty(branch)
        prepend!(ds1, rootcells)
        return true
    end
    # Merge branch into trunk
    for (i, cell) in enumerate(ds1)
        (alike, oppo) = cell
        if monoc(alike, branch)
            if isempty(oppo) || monoc(oppo, branch)
                verbose && println("Next cell")
                continue
            else
                append!(alike, branch)
                verbose && println("Merge alike")
                prepend!(ds1, rootcells)
                return true
            end
        else
            if isempty(oppo)
                # Collect higher trunk cells
                for (halike, hoppo) in ds1[i:end]
                    if !isempty(hoppo)
                        verbose && println("Not planer: dichromatic branch found")
                        return false
                    else
                        append!(oppo, halike)
                    end
                end
                newds = ds1[1:i-1]
                push!(newds, [branch, oppo])
                empty!(ds1)
                append!(ds1, newds)
                verbose && println("New opposite")
                prepend!(ds1, rootcells)
                return true
            elseif monoc(oppo, branch)
                append!(oppo, branch)
                verbose && println("Merge opposite")
                prepend!(ds1, rootcells)
                return true
            else
                verbose && println("Not planar: trichromatic")
                return false
            end
        end
    end
    verbose && println("Next cell")
    push!(ds1, Vector{Int}[branch, []])
    prepend!(ds1, rootcells)
    return true
end


function remove!(ds::DataStructure, edges)
    newds = DataStructure()
    for cell in ds
        (alike, oppo) = cell
        setdiff!(alike, edges)
        setdiff!(oppo, edges)
        if !isempty(alike)
            push!(newds, cell)
        end
    end
    empty!(ds)
    append!(ds, newds)
end


function planaritytest(graph::OrderedGraph; verbose=false)
    # Do DFS to determine treeedge, cotree, loworder, source
    state = DFSState(plaingraph(graph))
    dfs!(state)
    verbose && println(state)

    # L-R check
    ds = Dict(c => Vector{Vector{Int}}[[[c], []]] for c in keys(state.cotree))
    while !isempty(state.treeedge)
        e = pop!(state.treeedge)
        ds[e] = DataStructure()
        for i in state.loworder[e]
            verbose && println("stem: $(i), ds: $(ds[i])")
            if !merge!(ds[e], ds[i], state.cotree, verbose=verbose)
                empty!(state.treeedge)
                return false
            end
        end
        verbose && println("inedge: $(e), ds: $(ds[e])")
        finished = filter(
            x -> state.cotree[x] == state.source[e], keys(state.cotree))
        verbose && println("finished: $(finished)\n")
        remove!(ds[e], finished)
    end
    return true
end


function outerplanaritytest(graph::UndirectedGraph; verbose=false)
    # Add a node connects all other nodes
    newg = plaingraph(graph)
    newnode = addnode!(newg)
    for n in nodeset(graph)
        addedge!(newg, newnode, n)
    end
    return planaritytest(newg, verbose=verbose)
end


"""
    isplanar(graph::UndirectedGraph) -> Bool

Return whether the graph is planar.

# Reference

1. de Fraysseix, H., & Ossona de Mendez, P. (2012). Trémaux trees and planarity.
   European Journal of Combinatorics, 33(3), 279–293.
   https://doi.org/10.1016/j.ejc.2011.09.012
"""
@cache function isplanar(graph::UndirectedGraph)
    if isdefined(graph, :cache) && haskey(graph.cache, :isouterplanar)
        graph.cache[:isouterplanar] && return true
    end
    edgecount(graph) < 9 && return true
    edgecount(graph) > nodecount(graph) * 3 - 6 && return false
    return planaritytest(graph, verbose=verbose)
end


"""
    isouterplanar(graph::UndirectedGraph) -> Bool

Return whether the graph is outerplanar. The outerplanarity test is based on
a planarity test (see [`isplanar`](@ref)).
"""
@cache function isouterplanar(graph::UndirectedGraph)
    if isdefined(graph, :cache) && haskey(graph.cache, :isplanar)
        graph.cache[:isplanar] || return false
    end
    edgecount(graph) < 6 && return true
    edgecount(graph) > nodecount(graph) * 2 - 6 && return false
    return outerplanaritytest(graph, verbose=verbose)
end
