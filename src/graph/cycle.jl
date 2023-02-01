#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    mincyclebasis, mincycle_edges


function cotree_edges(g::SimpleGraph{T}) where T
    root = vertices(g)[1]
    stack = [root]
    visited = Set([root])
    stedges = Set{Edge{T}}()  # spanning tree edges
    while !isempty(stack)
        i = pop!(stack)
        for nbr in neighbors(g, i)
            if !(nbr in visited)
                push!(stedges, undirectededge(T, i, nbr))
                push!(visited, nbr)
                push!(stack, nbr)
            end
        end
    end
    return setdiff(Set(edges(g)), stedges)
end


function remap_edges(g::SimpleGraph{T}, func::Function) where T
    newadjdict = Dict{T,Vector{T}}()
    for (i, adjlist) in enumerate(g.fadjlist)
        for adj in adjlist
            u = func(i)
            haskey(newadjdict, u) || (newadjdict[u] = T[])
            push!(newadjdict[u], func(adj))
        end
    end
    nvx = maximum(keys(newadjdict))
    newadjlist = [T[] for _ in 1:nvx]
    for (k, v) in newadjdict
        push!(newadjlist[k], v...)
    end
    return SimpleGraph{T}(g.ne, newadjlist)
end


function disjointunion(g::SimpleGraph, h::SimpleGraph)
    h_ = remap_edges(h, i -> i + nv(g))
    hvmap = Dict(i + nv(g) => i for i in vertices(h))
    return union(g, h_), hvmap
end


function noweight_shortestpath(g::SimpleGraph{T}, u::T, v::T) where T
    # BFS based
    # u == v && return T[]  # impossible in findmincycle
    queue = [u]
    pred = Dict{Int,Int}(u => 0)
    while !isempty(queue)
        i = popfirst!(queue)
        for nbr in neighbors(g, i)
            if !haskey(pred, nbr)
                pred[nbr] = i
                push!(queue, nbr)
                if nbr == v  # target reached
                    empty!(queue)
                    break
                end
            end
        end
    end
    # retrieve path from the tree
    path = T[v]
    p = v
    while p != u
        p = pred[p]
        pushfirst!(path, p)
    end
    return path
end


function findmincycle(g::SimpleGraph, S::Set)
    subg, vmap = induced_subgraph(g, collect(setdiff(Set(edges(g)), S)))
    rev = Dict(v => i for (i, v) in enumerate(vmap))
    U, hvmap = disjointunion(subg, subg)
    hrev = Dict(v => k for (k, v) in hvmap)
    for s in S
        add_edge!(U, rev[src(s)], hrev[rev[dst(s)]])
        add_edge!(U, hrev[rev[src(s)]], rev[dst(s)])
    end
    paths = []
    for n in vertices(g)
        n1 = rev[n]
        n2 = hrev[rev[n]]
        sp = noweight_shortestpath(U, n1, n2)
        pop!(sp)
        push!(paths, sp)
    end
    return sortstablemin(paths, by=length)
end


"""
    mincyclebasis(graph::UndirectedGraph) -> Vector{Vector{Int}}

Returns minimum cycle basis represented as an array of edge sequence that make up a cycle.

# Reference

1. de Pina, J.C. Applications of shortest path methods. PhD thesis, University of Amsterdam,
   Nether-lands (1995)
1. Kavitha, T., Mehlhorn, K., Michail, D. & Paluch, K. E. An ``\\tilde{O}(m^{2}n)`` Algorithm for
   Minimum Cycle Basis of Graphs. Algorithmica 52, 333–349 (2008).
"""
function mincyclebasis(g::SimpleGraph{T}) where T
    cycles = Vector{T}[]
    for conn in connected_components(g)
        subg, vmap = induced_subgraph(g, conn)
        S = [Set([e]) for e in cotree_edges(subg)]
        N = length(S)  # N: circuit rank
        for k in 1:N
            p = findmincycle(subg, S[k])
            push!(cycles, p)  # TODO: apply vmap
            minedges = [undirectededge(T, p[i], p[i + 1]) for i in 1:(length(p) - 1)]
            for i in (k + 1):N
                if length(intersect(S[i], minedges)) % 2 == 1
                    S[i] = symdiff(S[i], S[k])
                end
            end
        end
    end
    return cycles
end


function mincycle_edges(g::SimpleGraph{T}) where T
    cycles = Vector{Edge{T}}[]
    for p in mincyclebasis(g)
        minedges = [undirectededge(T, p[i], p[i + 1]) for i in 1:(length(p) - 1)]
        push!(cycles, minedges)
    end
    return cycles
end