#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    shortestpath


function shortestpath(graph::AbstractUGraph, u, v)
    if u == v
        return
    end
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
