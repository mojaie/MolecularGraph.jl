#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    maxclique,
    maximalcliques


mutable struct FindCliqueState
    graph::UDGraph
    expire::UInt64
    Q::Set{Int}

    function FindCliqueState(g, timeout)
        exp = time_ns() + timeout * 1e9
        new(g, exp, Set())
    end
end


"""
    maxclique(graph::UDGraph; timeout::Int=3600) -> Set{Int}

Compute maximum clique of the graph. For details, see [`maximalcliques`](@ref).
"""
function maxclique(graph::UDGraph; timeout::Int=3600)
    res = maximalcliques(graph, timeout=timeout)
    # TODO: better way like python's max(iter, key=cmp)
    arr = collect(res)
    size = map(length, arr)
    return arr[argmax(size)]
end


"""
    maximalcliques(graph::UDGraph; timeout::Int=3600)

Return `Channel` which generates maximal cliques of the graph. Each cliques are
represented as a `Set` of member nodes.

# Reference
Tomita, E., Tanaka, A., & Takahashi, H. (2006). The worst-case time complexity
for generating all maximal cliques and computational experiments. Theoretical
Computer Science, 363(1), 28–42. https://doi.org/10.1016/J.TCS.2006.06.015
"""
function maximalcliques(graph::UDGraph; timeout::Int=3600)
    state = FindCliqueState(graph, timeout)
    return Channel(
        c::Channel -> expand!(state, nodekeys(graph), nodekeys(graph), c),
        ctype=Set{Int})
end


function expand!(state::FindCliqueState, subg, cand, channel)
    if isempty(subg)
        # Report max clique
        put!(channel, copy(state.Q))
        return
    elseif time_ns() > state.expire
        return
    end

    # Pivot
    # TODO: better way like python's max(iter, key=cmp)
    arr = collect(subg)
    deg = map(arr) do n
        length(intersect(cand, neighborkeys(state.graph, n)))
    end
    u = arr[argmax(deg)]

    unbrs = neighborkeys(state.graph, u)
    while !isempty(setdiff(cand, unbrs))
        q = pop!(setdiff(cand, unbrs))
        push!(state.Q, q)
        qnbrs = neighborkeys(state.graph, q)
        subgq = intersect(subg, qnbrs)
        candq = intersect(cand, qnbrs)
        expand!(state, subgq, candq, channel)
        pop!(cand, q, nothing)  # setdiff!(cand, Set([q]))
        # Revert
        pop!(state.Q, q)
    end
    return
end
