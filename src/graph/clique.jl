#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    maximumclique


mutable struct FindCliqueState{T<:UndirectedGraph}
    graph::T
    targetsize::Union{Int,Nothing} # Int

    adjacencies::Dict{Int,Set{Int}}
    expire::Union{UInt64,Nothing} # UInt64, nanoseconds
    Q::Set{Int}

    cliques::Vector{Set{Int}}
    status::Symbol

    function FindCliqueState{T}(graph; timeout=nothing, targetsize=nothing,
                kwargs...) where {T<:UndirectedGraph}
        if timeout !== nothing
            expire = (time_ns() + timeout * 1_000_000_000)::UInt64
        else
            expire = nothing
        end
        # fast adjacency access
        adj = Dict{Int,Set{Int}}()
        for n in nodeset(graph)
            adj[n] = adjacencies(graph, n)
        end
        new(graph, targetsize, adj, expire, Set(), [], :ready)
    end
end


function expand!(state::FindCliqueState, subg, cand)
    (state.status == :timedout || state.status == :targetreached) && return
    if isempty(subg)
        # Report max clique
        push!(state.cliques, copy(state.Q))
        return
    elseif state.expire !== nothing && time_ns() > state.expire
        state.status = :timedout
        return
    elseif state.targetsize !== nothing && length(state.Q) >= state.targetsize
        state.status = :targetreached
        push!(state.cliques, copy(state.Q))
        return
    end
    candnbrcnt(n) = length(intersect(cand, state.adjacencies[n]))
    pivot = sortstablemax(subg, by=candnbrcnt)

    for q in setdiff(cand, state.adjacencies[pivot])
        push!(state.Q, q)
        qnbrs = state.adjacencies[q]
        subgq = intersect(subg, qnbrs)
        candq = intersect(cand, qnbrs)
        expand!(state, subgq, candq)
        pop!(state.Q, q) # Revert
    end
    return
end


function expand!(
        state::FindCliqueState{ModularProduct}, subg, cand, qual)
    # c-clique version
    (state.status == :timedout || state.status == :targetreached) && return
    if isempty(subg)
        # Report max clique
        push!(state.cliques, copy(state.Q))
        return
    elseif state.expire !== nothing && time_ns() > state.expire
        state.status = :timedout
        return
    elseif state.targetsize !== nothing && length(state.Q) >= state.targetsize
        state.status = :targetreached
        push!(state.cliques, copy(state.Q))
        return
    end
    candnbrcnt(n) = length(intersect(cand, state.adjacencies[n]))
    pivot = sortstablemax(subg, by=candnbrcnt)

    for q in setdiff(cand, state.adjacencies[pivot])
        push!(state.Q, q)
        qconnbrs = Set{Int}()
        qdisnbrs = Set{Int}()
        for (inc, adj) in neighbors(state.graph, q)
            if edgeattr(state.graph, inc).hasedge
                push!(qconnbrs, adj)
            else
                push!(qdisnbrs, adj)
            end
        end
        qnbrs = state.adjacencies[q]
        subgq = union(intersect(subg, qnbrs), intersect(qual, qconnbrs))
        candq = union(intersect(cand, qnbrs), intersect(qual, qconnbrs))
        qualq = intersect(qual, qdisnbrs)
        expand!(state, subgq, candq, qualq)
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
    expand!(state, nodeset(graph), nodeset(graph))
    return sortstablemax(collect(state.cliques), by=length, init=[])
end
