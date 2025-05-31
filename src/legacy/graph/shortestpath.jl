#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    distance, reversedistance,
    reachablenodes, reversereachablenodes,
    isreachable, isreversereachable,
    eccentricity, reverseeccentricity, diameter,
    distancematrix,
    shortestpathnodes, reverseshortestpathnodes,
    shortestpathedges, reverseshortestpathedges,
    longestshortestpathnodes


function bfstree(adjfunc, graph, root)
    queue = [root]
    pred = Dict{Int,Int}(root => 0)
    while !isempty(queue)
        i = popfirst!(queue)
        for n in adjfunc(graph, i)
            if !haskey(pred, n)
                pred[n] = i
                push!(queue, n)
            end
        end
    end
    delete!(pred, root)
    return pred
end


function bfstree_edges(nbrfunc, graph, root)
    queue = [root]
    edges = Dict{Int,Int}(root => 0)
    while !isempty(queue)
        i = popfirst!(queue)
        for (inc, adj) in nbrfunc(graph, i)
            if !haskey(edges, adj)
                edges[adj] = inc
                push!(queue, adj)
            end
        end
    end
    delete!(edges, root)
    return edges
end


function bfsdepth(adjfunc, graph, root)
    queue = [root]
    depth = Dict{Int,Int}(root => 0)
    while !isempty(queue)
        i = popfirst!(queue)
        for n in adjfunc(graph, i)
            if !haskey(depth, n)
                depth[n] = depth[i] + 1
                push!(queue, n)
            end
        end
    end
    delete!(depth, root)
    return depth
end


"""
    distance(graph::AbstractGraph, source::Int, target::Int) -> Int
    reversedistance(graph::DirectedGraph, source::Int, target::Int) -> Int

Compute the distance (shortest path length) from `source` to `target`.
If the nodes are not reachable each other, the value will be `nothing`.
"""
function distance(adjfunc::Function, graph, source, target)
    source == target && return
    dps = bfsdepth(adjfunc, graph, source)
    haskey(dps, target) || return
    return dps[target]
end
distance(graph::UndirectedGraph, s, t) = distance(adjacencies, graph, s, t)
distance(graph::DirectedGraph, s, t) = distance(successors, graph, s, t)
reversedistance(graph::DirectedGraph, s, t) = distance(predecessors, graph, s, t)


"""
    reachablenodes(graph::AbstractGraph, node::Int) -> Set{Int}
    reversereachablenodes(graph::DirectedGraph, node::Int) -> Set{Int}

Return the set of reachable nodes from `node`.
"""
reachablenodes(graph::UndirectedGraph, node) = Set(keys(bfsdepth(adjacencies, graph, node)))
reachablenodes(graph::DirectedGraph, node) = Set(keys(bfsdepth(successors, graph, node)))
reversereachablenodes(graph::DirectedGraph, node) = Set(keys(bfsdepth(predecessors, graph, node)))


"""
    reachablenodes(graph::AbstractGraph, u::Int, v::Int) -> Bool
    reversereachablenodes(graph::DirectedGraph, u::Int, v::Int) -> Bool

Return whether the node `v` is reachable from `u`.
"""
isreachable(graph::UndirectedGraph, u, v) = distance(adjacencies, graph, u, v) !== nothing
isreachable(graph::DirectedGraph, u, v) = distance(successors, graph, u, v) !== nothing
isreversereachable(graph::DirectedGraph, u, v) = distance(predecessors, graph, u, v) !== nothing


"""
    eccentricity(graph::UndirectedGraph, v::Int) -> Int

Compute the eccentricity of the graph (the largest distance between `v` and
any other nodes).
"""
function eccentricity(adjfunc::Function, graph, v)
    dists = values(bfsdepth(adjfunc, graph, v))
    isempty(dists) && return
    return maximum(dists)
end
eccentricity(graph::UndirectedGraph, v) = eccentricity(adjacencies, graph, v)
eccentricity(graph::DirectedGraph, v) = eccentricity(successors, graph, v)
reverseeccentricity(graph::DirectedGraph, v) = eccentricity(predecessors, graph, v)


"""
    diameter(graph::AbstractGraph) -> Int

Compute the diameter of the graph (the largest eccentricity of any nodes).
"""
function diameter(graph::AbstractGraph)
    nodecount(graph) == 0 && return
    eccs = Int[]
    for n in nodeset(graph)
        ecc = eccentricity(graph, n)
        ecc === nothing || push!(eccs, ecc)
    end
    return maximum(eccs)
end


"""
    distancematrix(graph::OrderedGraph) -> Matrix{Float64}

Generate the distance matrix of the graph.

Note that the type of the generated matrix will be `Float64`. If the nodes are
not reachable each other, the distance value will be `Inf`.
"""
function distancematrix(adjfunc::Function, graph)
    ncnt = nodecount(graph)
    distmat = fill(Inf, (ncnt, ncnt))
    for s in 1:ncnt
        for (t, d) in bfsdepth(adjfunc, graph, s)
            distmat[s, t] = d
        end
    end
    return distmat
end

distancematrix(graph::OrderedGraph) = distancematrix(adjacencies, graph)
distancematrix(graph::OrderedDiGraph) = distancematrix(successors, graph)


"""
    shortestpathnodes(graph::UndirectedGraph, u::Int, v::Int) -> Vector{Int}

Compute the shortest path between `u` and `v` as a vector of the nodes that
forms the path. Return nothing if not reachable.
"""
function shortestpathnodes(func::Function, graph, source, target)
    source == target && return
    tree = bfstree(func, graph, source)
    path = Int[target]
    p = target
    while p != source
        p = tree[p]
        pushfirst!(path, p)
    end
    return path
end

shortestpathnodes(graph::UndirectedGraph, u, v) = shortestpathnodes(adjacencies, graph, u, v)
shortestpathnodes(graph::DirectedGraph, u, v) = shortestpathnodes(successors, graph, u, v)
reverseshortestpathnodes(graph::DirectedGraph, u, v) = shortestpathnodes(predecessors, graph, u, v)


"""
    shortestpathedges(graph::UndirectedGraph, u::Int, v::Int) -> Vector{Int}

Compute the shortest path between `u` and `v` as a vector of the edges that
forms the path. Return nothing if not reachable.
"""
function shortestpathedges(func::Function, graph, source, target)
    source == target && return
    tree = bfstree_edges(func, graph, source)
    path = Int[]
    p = target
    while p != source
        inc = tree[p]
        p = neighbors(graph, p)[inc]
        pushfirst!(path, inc)
    end
    return path
end

shortestpathedges(graph::UndirectedGraph, u, v) = shortestpathedges(neighbors, graph, u, v)
shortestpathedges(graph::DirectedGraph, u, v) = shortestpathedges(outneighbors, graph, u, v)
reverseshortestpathedges(graph::DirectedGraph, u, v) = shortestpathedges(inneighbors, graph, u, v)


"""
    longestshortestpathnodes(graph::UndirectedGraph) -> Vector{Int}

Compute the longest shortest path in the graph (a path between two arbitrary
peripheral nodes) as a vector of nodes that starts with one of the
peripheral node and ends with the other side.
"""
function longestshortestpathnodes(graph::UndirectedGraph)
    peribk = Dict{Int,Tuple{Int,Int}}()
    for u in nodeset(graph)
        dps = bfsdepth(adjacencies, graph, u)
        isempty(dps) && continue
        # bucket sort
        bk = Dict{Int,Int}()
        for (k, v) in dps
            bk[v] = k
        end
        maxlen = maximum(keys(bk))
        v = bk[maxlen]
        peribk[maxlen] = (u, v)
    end
    isempty(peribk) && return
    peripheral = peribk[maximum(keys(peribk))]
    return shortestpathnodes(graph, peripheral...)
end
