#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    bfs_path, noweight_shortestpath, noweight_all_distances


function bfs_path(parents, t)
    path = [t]
    x = t
    while parents[x] != x
        x = parents[x]
        pushfirst!(path, x)
    end
    return path
end


noweight_shortestpath(g::SimpleGraph, s, t) = bfs_path(bfs_parents(g, s), t)


function noweight_all_distances(g::SimpleGraph{T}, s, t) where T
    """ DFS based algorithm to enumerate all distances of
    the paths between two vertices."""
    stack = [s]
    parent = Dict(s => s)
    cycle = Dict{T,Vector{T}}(i => T[] for i in vertices(g))
    tpaths = []
    while !isempty(stack)
        i = pop!(stack)
        for nbr in neighbors(g, i)
            nbr == parent[i] && continue
            if nbr == t  # target found
                dp = bfs_path(parent, i)
                for v in dp
                    push!(cycle[v], length(tpaths) + 1)
                end
                push!(tpaths, dp)
            elseif nbr in keys(parent)  # cycle found
                dp = bfs_path(parent, i)
                for n in copy(cycle[nbr])
                    tp = tpaths[n]
                    c = findfirst(x -> x == nbr, tp)
                    merged = vcat(dp, tp[c:end])
                    for v in merged
                        push!(cycle[v], length(tpaths) + 1)
                    end
                    push!(tpaths, merged)
                end
            else  # new vertex
                parent[nbr] = i
                push!(stack, nbr)
            end
        end
    end
    return sort(collect(Set(map(length, tpaths))))
end