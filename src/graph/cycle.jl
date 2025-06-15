#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

function cotree_edges(g::SimpleGraph{T}) where T
    root = vertices(g)[1]
    stack = [root]
    visited = Set([root])
    stedges = Set{Edge{T}}()  # spanning tree edges
    while !isempty(stack)
        i = pop!(stack)
        for nbr in neighbors(g, i)
            if !(nbr in visited)
                push!(stedges, u_edge(T, i, nbr))
                push!(visited, nbr)
                push!(stack, nbr)
            end
        end
    end
    return setdiff(Set(edges(g)), stedges)
end


function findmincycle(g::SimpleGraph, S::Set)
    subg, vmap = induced_subgraph(g, collect(setdiff(Set(edges(g)), S)))
    rev = Dict(v => i for (i, v) in enumerate(vmap))
    U = disjoint_union(subg, subg)
    hvmap = Dict(i + nv(g) => i for i in vertices(subg))
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
    # resume original vertices and return the cycle
    return [vmap[get(hvmap, s, s)] for s in sortstablemin(paths, by=length)]
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
        @debug "init" S
        N = length(S)  # N: circuit rank
        for k in 1:N
            @debug "S" S[k]
            p = findmincycle(subg, S[k])
            @debug "mincycle" vmap[p]
            push!(cycles, vmap[p])
            minedges = [u_edge(T, p[i], p[i + 1]) for i in 1:(length(p) - 1)]
            push!(minedges, u_edge(T, p[1], p[length(p)]))
            for i in (k + 1):N
                if length(intersect(S[i], minedges)) % 2 == 1
                    @debug "symdiff" i symdiff(S[i], S[k])
                    S[i] = symdiff(S[i], S[k])
                end
            end
        end
    end
    return cycles
end


function edgemincyclebasis(g::SimpleGraph{T}) where T
    cycles = Vector{Edge{T}}[]
    for p in mincyclebasis(g)
        minedges = [u_edge(T, p[i], p[i + 1]) for i in 1:(length(p) - 1)]
        push!(minedges, u_edge(T, p[1], p[end]))
        push!(cycles, minedges)
    end
    return cycles
end
