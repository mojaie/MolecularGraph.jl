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
    u::Set{Int}
    v::Set{Int}
    adj::Dict{Int,Set{Int}}
    dist::Dict{Union{Int,Nothing},Union{Int,Nothing}}
    utov::Dict{Int,Union{Int,Nothing}}
    vtou::Dict{Int,Union{Int,Nothing}}

    function MaxCardMatchState(u, v, adj)
        new(Set(u), Set(v), adj, Dict(),
            Dict(i => nothing for i in u), Dict(i => nothing for i in v))
    end
end


function MaxCardMatchState(u, v, efunc::Function)
    # Generate bipartite relationship
    adj = Dict{Int,Set{Int}}()
    for i in u
        adj[i] = Set{Int}()
        for j in v
            if efunc(i, j)
                push!(adj[u], j)
            end
        end
    end
    MaxCardMatchState(Set(u), Set(v), adj)
end


function bfs!(state::MaxCardMatchState)
    queue = []
    for i in state.u
        if state.utov[i] === nothing
            state.dist[i] = 0
            push!(queue, i)
        else
            state.dist[i] = nothing
        end
    end
    state.dist[nothing] = nothing
    while !isempty(queue)
        i = popfirst!(queue)
        ndist = state.dist[nothing]
        if ndist === nothing || state.dist[i] < ndist
            for j in state.adj[i]
                m = state.vtou[j]
                if state.dist[m] === nothing
                    state.dist[m] = state.dist[i] + 1
                    push!(queue, m)
                end
            end
        end
    end
    return state.dist[nothing] !== nothing
end


function dfs!(state::MaxCardMatchState, i)
    if i === nothing
        return true
    end
    for j in state.adj[i]
        m = state.vtou[j]
        if state.dist[m] == state.dist[i] + 1
            if dfs!(state, m)
                state.vtou[j] = i
                state.utov[i] = j
                return true
            end
        end
    end
    state.dist[i] = nothing
    return false
end


function maxcardmatch!(state::MaxCardMatchState)
    # Hopcroft-Karp method
    while bfs!(state)
        for i in state.u
            if state.utov[i] === nothing
                dfs!(state, i)
            end
        end
    end
end


function maxcardmap(u, v, adj)
    # Maximum cardinality mapping
    state = MaxCardMatchState(u, v, adj)
    maxcardmatch!(state)
    mapping = Dict{Int,Int}()
    for (i, j) in state.utov
        if i !== nothing && j !== nothing
            mapping[i] = j
        end
    end
    return mapping
end

maxcard(u, v, adj) = length(maxcardmap(u, v, adj))