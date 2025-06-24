#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

function bfs_path(parents, t::T) where T<:Integer
    path = [t]
    x = t
    while parents[x] != x
        x = parents[x]
        pushfirst!(path, x)
    end
    return path
end


function unweighted_shortestpath(g::SimpleGraph{T}, src::T, dst::T) where T <: Integer
    # Specialized BFS implementation for molecular graphs.
    # Reference: https://github.com/JuliaGraphs/Graphs.jl/blob/master/src/traversals/bfs.jl
    n = nv(g)
    dst > n && return T[]  # Invalid target
    parents = zeros(T, n)
    queue = [src]
    parents[src] = src
    while !isempty(queue)
        v = popfirst!(queue)
        for nbr in neighbors(g, v)
            if parents[nbr] == 0
                push!(queue, nbr)
                parents[nbr] = v
            end
            nbr == dst && return bfs_path(parents, nbr)  # target reached
        end
    end
    return T[]  # Unreachable target
end


function noweight_all_distances(g::SimpleGraph{T}, s::T, t::T) where T<:Integer
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