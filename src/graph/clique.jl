#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    maximumclique,
    maximalcliques


mutable struct FindCliqueState{T<:UndirectedGraph}
    # Input
    graph::T
    # Optional
    timeout # Int
    c_clique_constraint # Function
    # State
    Q::Set{Int}
    expire # nanoseconds, UInt64
    status::Symbol

    function FindCliqueState{T}(graph) where {T<:UndirectedGraph}
        state = new()
        state.graph = graph
        state.Q = Set()
        state.status = :Ready
        return state
    end
end


"""
    maximumclique(graph::UndirectedGraph; kwargs...) -> Set{Int}

Compute maximum clique of the graph. For details, see [`maximalcliques`](@ref).
"""
maximumclique(
    graph::UndirectedGraph; kwargs...
) = sortstablemax(
    collect(maximalcliques(graph; kwargs...)), by=length, init=[])


"""
    maximalcliques(graph::UndirectedGraph; kwargs...)

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
function maximalcliques(graph::T; kwargs...) where {T<:UndirectedGraph}
    state = FindCliqueState{T}(graph)
    if haskey(kwargs, :timeout)
        state.timeout = kwargs[:timeout]::Int
        state.expire = (time_ns() + state.timeout * 1_000_000_000)::UInt64
    end
    if haskey(kwargs, :c_clique_constraint)
        state.c_clique_constraint = kwargs[:c_clique_constraint]::Function
        f = c::Channel -> expand!(
            state, nodeset(graph), nodeset(graph), nodeset(graph), c)
    else
        f = c::Channel -> expand!(state, nodeset(graph), nodeset(graph), c)
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
    candnbrcnt(n) = length(intersect(cand, adjacencies(state.graph, n)))
    pivot = sortstablemax(subg, by=candnbrcnt)

    for q in setdiff(cand, adjacencies(state.graph, pivot))
        push!(state.Q, q)
        qnbrs = adjacencies(state.graph, q)
        subgq = intersect(subg, qnbrs)
        candq = intersect(cand, qnbrs)
        expand!(state, subgq, candq, channel)
        pop!(state.Q, q) # Revert
    end
    return
end


function expand!(
        state::FindCliqueState{ModularProduct}, subg, cand, qual, channel)
    # c-clique version
    if isempty(subg)
        # Report max clique
        put!(channel, copy(state.Q))
        return
    elseif isdefined(state, :timeout) && time_ns() > state.expire
        return
    end
    candnbrcnt(n) = length(intersect(cand, adjacencies(state.graph, n)))
    pivot = sortstablemax(subg, by=candnbrcnt)

    for q in setdiff(cand, adjacencies(state.graph, pivot))
        push!(state.Q, q)
        qconnbrs = Set{Int}()
        qdisnbrs = Set{Int}()
        for (inc, adj) in neighbors(state.graph, q)
            if getedgeattr(state.graph, inc).hasedge
                push!(qconnbrs, adj)
            else
                push!(qdisnbrs, adj)
            end
        end
        qnbrs = adjacencies(state.graph, q)
        subgq = union(intersect(subg, qnbrs), intersect(qual, qconnbrs))
        candq = union(intersect(cand, qnbrs), intersect(qual, qconnbrs))
        qualq = intersect(qual, qdisnbrs)
        expand!(state, subgq, candq, qualq, channel)
        pop!(state.Q, q) # Revert
    end
    return
end
