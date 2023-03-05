#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    MaxCardMatchState,
    maxcardmap,
    maxcard,
    maxcardmatch!


mutable struct MaxCardMatchState
    U::Set{Int}
    V::Set{Int}
    adj::Dict{Int,Set{Int}}
    dist::Dict{Union{Int,Nothing},Union{Int,Nothing}}
    utov::Dict{Int,Union{Int,Nothing}}
    vtou::Dict{Int,Union{Int,Nothing}}

    function MaxCardMatchState(U, V, adj::Dict)
        new(Set(U), Set(V), adj, Dict(),
            Dict(u => nothing for u in U), Dict(v => nothing for v in V))
    end
end


function MaxCardMatchState(U, V, efunc::Function)
    # Generate bipartite relationship
    adj = Dict{Int,Set{Int}}()
    for u in U
        adj[u] = Set{Int}()
        for v in V
            if efunc(u, v)
                push!(adj[u], v)
            end
        end
    end
    MaxCardMatchState(Set(U), Set(V), adj)
end


function maxcardmap(U, V, adj)
    # Maximum cardinality mapping
    state = MaxCardMatchState(U, V, adj)
    maxcardmatch!(state)
    mapping = Dict{Int,Int}()
    for (u, v) in state.utov
        if u !== nothing && v !== nothing
            mapping[u] = v
        end
    end
    return mapping
end

maxcard(U, V, adj) = length(maxcardmap(U, V, adj))


function maxcardmatch!(state::MaxCardMatchState)
    # Hopcroft-Karp method
    while bfs!(state)
        for u in state.U
            if state.utov[u] === nothing
                dfs!(state, u)
            end
        end
    end
end


function bfs!(state::MaxCardMatchState)
    queue = []
    for u in state.U
        if state.utov[u] === nothing
            state.dist[u] = 0
            push!(queue, u)
        else
            state.dist[u] = nothing
        end
    end
    state.dist[nothing] = nothing
    while !isempty(queue)
        u = popfirst!(queue)
        ndist = state.dist[nothing]
        if ndist === nothing || state.dist[u] < ndist
            for v in state.adj[u]
                m = state.vtou[v]
                if state.dist[m] === nothing
                    state.dist[m] = state.dist[u] + 1
                    push!(queue, m)
                end
            end
        end
    end
    return state.dist[nothing] !== nothing
end


function dfs!(state::MaxCardMatchState, u)
    if u === nothing
        return true
    end
    for v in state.adj[u]
        m = state.vtou[v]
        if state.dist[m] == state.dist[u] + 1
            if dfs!(state, m)
                state.vtou[v] = u
                state.utov[u] = v
                return true
            end
        end
    end
    state.dist[u] = nothing
    return false
end
