#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    maximumclique


mutable struct FindCliqueState{T<:UndirectedGraph}
    graph::T
    timeout # Int
    targetsize # Int
    c_clique_constraint # Function

    expire # nanoseconds, UInt64

    Q::Set{Int}
    channel::Channel
    status::Symbol

    function FindCliqueState{T}(graph; kwargs...) where {T<:UndirectedGraph}
        state = new()
        state.graph = graph
        state.Q = Set()
        state.status = :ready
        if haskey(kwargs, :timeout)
            state.timeout = kwargs[:timeout]::Int
            state.expire = (time_ns() + state.timeout * 1_000_000_000)::UInt64
        end
        if haskey(kwargs, :targetsize)
            state.targetsize = kwargs[:targetsize]
        end
        if haskey(kwargs, :c_clique_constraint)
            state.c_clique = kwargs[:c_clique_constraint]::Function
            ch = c::Channel -> expand!(
                state, nodeset(graph), nodeset(graph), nodeset(graph), c)
        else
            ch = c::Channel -> expand!(
                state, nodeset(graph), nodeset(graph), c)
        end
        state.channel = Channel(ch, ctype=Set{Int})
        return state
    end
end


function expand!(state::FindCliqueState, subg, cand, channel)
    if isempty(subg)
        # Report max clique
        put!(channel, copy(state.Q))
        return
    elseif isdefined(state, :timeout) && time_ns() > state.expire
        state.status = :timedout
        return
    elseif isdefined(state, :targetsize) && length(state.Q) > state.targetsize
        state.status = :targetreached
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
        state.status = :timedout
        return
    elseif isdefined(state, :targetsize) && length(state.Q) > state.targetsize
        state.status = :targetreached
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


"""
    maximumclique(graph::UndirectedGraph; kwargs...) -> Set{Int}

Return a set of maximum clique nodes.

# Reference

1. Tomita, E., Tanaka, A., & Takahashi, H. (2006). The worst-case time
   complexity for generating all maximal cliques and computational experiments.
   Theoretical Computer Science, 363(1), 28–42.
   https://doi.org/10.1016/J.TCS.2006.06.015
1. Cazals, F., & Karande, C. (2005). An algorithm for reporting maximal
   c-cliques. Theoretical Computer Science, 349(3), 484–490.
   https://doi.org/10.1016/j.tcs.2005.09.038

"""
function maximumclique(graph::T; kwargs...) where {T<:UndirectedGraph}
    state = FindCliqueState{T}(graph; kwargs...)
    return sortstablemax(collect(state.channel), by=length, init=[])
end
