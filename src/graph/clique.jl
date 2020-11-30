#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    maximalcliques, maximumclique,
    maximalconncliques, maximumconnclique


mutable struct FindCliqueState{T<:UndirectedGraph}
    graph::T
    targetsize::Union{Int,Nothing} # Int

    adjacencies::Dict{Int,Set{Int}}
    expire::Union{UInt64,Nothing} # UInt64, nanoseconds
    Q::Vector{Int}

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
        new(graph, targetsize, adj, expire, [], [], :ongoing)
    end
end
FindCliqueState(graph::UndirectedGraph; kwargs...) = FindCliqueState{typeof(graph)}(graph; kwargs...)


mutable struct FindConnCliqueState{T<:ModularProduct}
    graph::T
    targetsize::Union{Int,Nothing} # Int

    adjacencies::Dict{Int,Set{Int}}
    connected::Dict{Int,Set{Int}}
    disconn::Dict{Int,Set{Int}}
    expire::Union{UInt64,Nothing} # UInt64, nanoseconds

    cliques::Vector{Set{Int}}
    status::Symbol

    function FindConnCliqueState{T}(graph; timeout=nothing, targetsize=nothing,
                kwargs...) where {T<:ModularProduct}
        if timeout !== nothing
            expire = (time_ns() + timeout * 1_000_000_000)::UInt64
        else
            expire = nothing
        end
        # fast adjacency access
        adj = Dict{Int,Set{Int}}()
        conn = Dict{Int,Set{Int}}()
        disconn = Dict{Int,Set{Int}}()
        for n in nodeset(graph)
            adj[n] = Set{Int}()
            conn[n] = Set{Int}()
            disconn[n] = Set{Int}()
            for (i, a) in neighbors(graph, n)
                push!(adj[n], a)
                if edgeattr(graph, i).hasedge
                    push!(conn[n], a)
                else
                    push!(disconn[n], a)
                end
            end
        end
        new(graph, targetsize, adj, conn, disconn, expire, [], :ongoing)
    end
end
FindConnCliqueState(graph::ModularProduct; kwargs...) = FindConnCliqueState{typeof(graph)}(graph; kwargs...)



function expand!(state::FindCliqueState, subg, cand)
    (state.status == :timedout || state.status == :targetreached) && return
    if isempty(subg)
        # Report max clique
        push!(state.cliques, Set(state.Q))
        return
    elseif state.expire !== nothing && time_ns() > state.expire
        state.status = :timedout
        return
    elseif state.targetsize !== nothing && length(state.Q) >= state.targetsize
        state.status = :targetreached
        push!(state.cliques, Set(state.Q))
        return
    end
    candnbrcnt(n) = length(intersect(cand, state.adjacencies[n]))
    pivot = sortstablemax(subg, by=candnbrcnt)
    copv = setdiff(cand, state.adjacencies[pivot])

    for q in copv
        push!(state.Q, q)
        qnbrs = state.adjacencies[q]
        subgq = intersect(subg, qnbrs)
        candq = intersect(cand, qnbrs)
        expand!(state, subgq, candq)
        pop!(cand, q)
        pop!(state.Q)
    end
    return
end

function expandconn!(state::FindConnCliqueState, R, P, Q, X, Y)
    (state.status == :timedout || state.status == :targetreached) && return
    if isempty(P) && isempty(X)
        # Report max clique
        push!(state.cliques, copy(R))
        return
    elseif state.expire !== nothing && time_ns() > state.expire
        state.status = :timedout
        return
    elseif state.targetsize !== nothing && length(R) >= state.targetsize
        state.status = :targetreached
        push!(state.cliques, copy(R))
        return
    end
    while !isempty(P)
        n = pop!(P)
        Rnew = union(R, [n])
        Qnew = intersect(Q, state.disconn[n])
        Pnew = union(
            intersect(P, state.adjacencies[n]),
            intersect(Q, state.connected[n]))
        Ynew = intersect(Y, state.disconn[n])
        Xnew = union(
            intersect(X, state.adjacencies[n]),
            intersect(Y, state.connected[n]))
        expandconn!(state, Rnew, Pnew, Qnew, Xnew, Ynew)
        push!(X, n)
    end
    return
end


"""
    maximalcliques(graph::UndirectedGraph; kwargs...
        ) -> Tuple{Vector{Set{Int}}, Symbol}

Return maximal cliques.

# Reference

1. Tomita, E., Tanaka, A., & Takahashi, H. (2006). The worst-case time
   complexity for generating all maximal cliques and computational experiments.
   Theoretical Computer Science, 363(1), 28–42.
   https://doi.org/10.1016/J.TCS.2006.06.015
1. Cazals, F., & Karande, C. (2008). A note on the problem of reporting maximal
   cliques. Theoretical Computer Science, 407(1–3), 564–568.
   https://doi.org/10.1016/j.tcs.2008.05.010
"""
function maximalcliques(graph::UndirectedGraph; kwargs...)
    state = FindCliqueState(graph; kwargs...)
    expand!(state, nodeset(graph), nodeset(graph))
    if state.status == :ongoing
        state.status = :done
    end
    return (state.cliques, state.status)
end


"""
    maximumclique(graph::UndirectedGraph; kwargs...) -> Tuple{Set{Int}, Symbol}

Return a maximum clique.

"""
function maximumclique(graph::UndirectedGraph; kwargs...)
    (cliques, status) = maximalcliques(graph; kwargs...)
    return (sortstablemax(cliques, by=length, init=[]), status)
end



"""
    maximalconncliques(graph::ModularProduct; kwargs...
        ) -> Tuple{Vector{Set{Int}}, Symbol}

Return maximal connected cliques.

# Reference

1. Cazals, F., & Karande, C. (2005). An algorithm for reporting maximal
   c-cliques. Theoretical Computer Science, 349(3), 484–490.
   https://doi.org/10.1016/j.tcs.2005.09.038

"""
function maximalconncliques(graph::ModularProduct; kwargs...)
    state = FindConnCliqueState(graph; kwargs...)
    nodes = nodeset(graph)
    done = Set{Int}()
    for n in nodes
        R = Set([n])
        P = intersect(setdiff(nodes, done), state.connected[n])
        Q = intersect(setdiff(nodes, done), state.disconn[n])
        X = intersect(state.connected[n], done)
        Y = intersect(state.disconn[n], done)
        expandconn!(state, R, P, Q, X, Y)
        push!(done, n)
    end
    if state.status == :ongoing
        state.status = :done
    end
    return (state.cliques, state.status)
end


"""
    maximumconnclique(graph::ModularProduct; kwargs...
        ) -> Tuple{Set{Int}, Symbol}

Return a maximum connected clique.

"""
function maximumconnclique(graph::ModularProduct; kwargs...)
    (cliques, status) = maximalconncliques(graph; kwargs...)
    return (sortstablemax(cliques, by=length, init=[]), status)
end
