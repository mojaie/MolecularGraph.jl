#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    shortestpath,
    distance, distancematrix,
    eccentricity, diameter,
    longestshortestpath,
    reachablenodes, pathlength, backtrack


"""
    shortestpath(graph::UDGraph, u, v) -> Vector{Int}

Compute the shortest path between `u` and `v` as a vector of the nodes that
starts with `u` and ends with`v`. Return nothing if u == v or not reachable.
"""
function shortestpath(graph::UDGraph, u, v)
    u == v && return
    queue = [u]
    pred = Dict{Int,Union{Int,Nothing}}(u => nothing)
    while !isempty(queue)
        i = popfirst!(queue)
        for nbr in neighborkeys(graph, i)
            if nbr == v
                # v found
                path = [v]
                p = i
                while p !== nothing
                    pushfirst!(path, p)
                    p = pred[p]
                end
                return path
            elseif !(nbr in keys(pred))
                # New nodes
                pred[nbr] = i
                push!(queue, nbr)
            end
        end
    end
    # return nothing if not reachable
end


"""
    distance(graph::UDGraph, root) -> Dict{Int,Union{Int,Nothing}}

Compute the distance from `root` to any other nodes. If the nodes are not
reachable each other, the value will be `nothing`.
"""
function distance(graph::UDGraph, root)
    queue = [root]
    dist = Dict{Int,Union{Int,Nothing}}(root => 0)
    while !isempty(queue)
        q = popfirst!(queue)
        for nbr in neighborkeys(graph, q)
            if !(nbr in keys(dist))
                dist[nbr] = dist[q] + 1
                push!(queue, nbr)
            end
        end
    end
    return dist
end


"""
    distancematrix(graph::UDGraph) -> Matrix{Float64}

Compute the distance among each other nodes.

Note that the type of the generated matrix will be `Float64`. If the nodes are
not reachable each other, the distance value will be `Inf`.
"""
function distancematrix(graph::UDGraph)
    ncnt = nodecount(graph)
    distmat = zeros((ncnt, ncnt))
    for s in nodekeys(graph) # TODO: should be ordered
        for (t, d) in distance(graph, s)
            distmat[s, t] = d === nothing ? Inf : d
        end
    end
    return distmat
end


"""
    eccentricity(graph::UDGraph, v) => Int

Compute the eccentricity of the graph (the largest distance between `v` and
any other nodes).
"""
eccentricity(graph::UDGraph, v) = maximum(values(distance(graph, v)))


"""
    diameter(graph::UDGraph) => Int

Compute the diameter of the graph (the largest eccentricity of any nodes).
"""
diameter(graph::UDGraph) = maximum(
    eccentricity(graph, n) for n in nodekeys(graph))


"""
    longestshortestpath(graph::UDGraph) -> Vector{Int}

Compute the longest shortest path in the graph (a path between two arbitrary
peripheral nodes) as a vector of nodes that starts with one of the
peripheral node and ends with the other side.
"""
function longestshortestpath(graph::UDGraph)
    otherside = Dict{Int,Int}()
    eccent = Dict{Int,Int}()
    for v in nodekeys(graph)
        d = distance(graph, v)
        farthest = sort(collect(keys(d)), by=k->d[k])[end]
        otherside[v] = farthest
        eccent[v] = d[farthest]
    end
    peripheral = sort(collect(keys(eccent)), by=k->d[k])[end]
    return shortestpath(graph, peripheral, otherside[peripheral])
end



# TODO: refactor DGraph methods


function reachablenodes(graph::DGraph, source)
    pred = shortestpath(graph, source)
    delete!(pred, source)
    return keys(pred)
end


function pathlength(graph::DGraph, source)
    res = Dict{Int,Int}()
    pred = shortestpath(graph, source)
    for p in keys(pred)
        path = backtrack(pred, source, p)
        res[p] = length(path) - 1
    end
    return res
end


function shortestpath(graph::DGraph, source)
    queue = [source]
    pred = Dict{Int,Union{Int,Nothing}}(source => nothing)
    while !isempty(queue)
        i = popfirst!(queue)
        for succ in succkeys(graph, i)
            if !(succ in keys(pred))
                # New nodes
                pred[succ] = i
                push!(queue, succ)
            end
        end
    end
    return pred
end


function backtrack(pred::Dict{Int,Union{Int,Nothing}}, source, target)
    path = []
    p = target
    while p !== nothing
        pushfirst!(path, p)
        if p == source
            break
        end
        p = pred[p]
    end
    return path
end
