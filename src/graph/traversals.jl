#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    bfs_path, noweight_shortestpath


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
