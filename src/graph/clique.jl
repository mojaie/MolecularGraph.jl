#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    maxclique,
    maximalcliques


mutable struct FindCliqueState
    # Input
    graph::UDGraph
    # Optional
    timeout # Int
    c_clique_constraint # Function
    # State
    Q::Set{Int}
    expire # nanoseconds, UInt64
    status::Symbol

    function FindCliqueState(graph)
        state = new()
        state.graph = graph
        state.Q = Set()
        state.status = :Ready
        return state
    end
end


"""
    maxclique(graph::UDGraph; timeout::Int=3600) -> Set{Int}

Compute maximum clique of the graph. For details, see [`maximalcliques`](@ref).
"""
function maxclique(graph::UDGraph; kwargs...)
    # TODO: better way like python's max(iter, key=cmp)
    maxclq = Set{Int}()
    for c in maximalcliques(graph; kwargs...)
        if length(c) > length(maxclq)
            maxclq = c
        end
    end
    return maxclq
end


"""
    maximalcliques(graph::UDGraph; kwargs...)

Return `Channel` which generates maximal cliques of the graph. Each cliques are
represented as a `Set` of member nodes.

# Reference

1. Tomita, E., Tanaka, A., & Takahashi, H. (2006). The worst-case time
complexity for generating all maximal cliques and computational experiments.
Theoretical Computer Science, 363(1), 28–42.
https://doi.org/10.1016/J.TCS.2006.06.015

1. Cazals, F., & Karande, C. (2005). An algorithm for reporting maximal
c-cliques. Theoretical Computer Science, 349(3), 484–490.
https://doi.org/10.1016/j.tcs.2005.09.038

"""
function maximalcliques(graph::UDGraph; kwargs...)
    state = FindCliqueState(graph)
    if haskey(kwargs, :timeout)
        state.timeout = kwargs[:timeout]::Int
        state.expire = (time_ns() + state.timeout * 1_000_000_000)::UInt64
    end
    if haskey(kwargs, :c_clique_constraint)
        state.c_clique_constraint = kwargs[:c_clique_constraint]::Function
        f = c::Channel -> expand!(
            state, nodekeys(graph), nodekeys(graph), nodekeys(graph), c)
    else
        f = c::Channel -> expand!(state, nodekeys(graph), nodekeys(graph), c)
    end
    return Channel(f, ctype=Set{Int})
end


function expand!(state::FindCliqueState, subg, cand, channel)
    if isempty(subg)
        # Report max clique
        put!(channel, copy(state.Q))
        return
    elseif isdefined(state, :timeout) && time_ns() > state.expire
        return
    end
    # Pivot
    # TODO: better way like python's max(iter, key=cmp)
    arr = collect(subg)
    deg = map(arr) do n
        length(intersect(cand, neighborkeys(state.graph, n)))
    end
    pv = arr[argmax(deg)]

    for q in setdiff(cand, neighborkeys(state.graph, pv))
        push!(state.Q, q)
        qnbrs = neighborkeys(state.graph, q)
        subgq = intersect(subg, qnbrs)
        candq = intersect(cand, qnbrs)
        expand!(state, subgq, candq, channel)
        pop!(state.Q, q) # Revert
    end
    return
end


function expand!(state::FindCliqueState, subg, cand, qual, channel)
    # c-clique version
    if isempty(subg)
        # Report max clique
        put!(channel, copy(state.Q))
        return
    elseif isdefined(state, :timeout) && time_ns() > state.expire
        return
    end
    # Pivot
    # TODO: better way like python's max(iter, key=cmp)
    arr = collect(subg)
    deg = map(arr) do n
        length(intersect(cand, neighborkeys(state.graph, n)))
    end
    pv = arr[argmax(deg)]

    for q in setdiff(cand, neighborkeys(state.graph, pv))
        push!(state.Q, q)
        qconnbrs = Set{Int}()
        qdisnbrs = Set{Int}()
        for (n, nbr) in neighbors(state.graph, q)
            if getedge(state.graph, nbr).hasedge
                push!(qconnbrs, n)
            else
                push!(qdisnbrs, n)
            end
        end
        qnbrs = neighborkeys(state.graph, q)
        subgq = union(intersect(subg, qnbrs), intersect(qual, qconnbrs))
        candq = union(intersect(cand, qnbrs), intersect(qual, qconnbrs))
        qualq = intersect(qual, qdisnbrs)
        expand!(state, subgq, candq, qualq, channel)
        pop!(state.Q, q) # Revert
    end
    return
end
