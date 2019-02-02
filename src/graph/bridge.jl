#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    bridges


mutable struct FindBridgeState{G<:UDGraph}
    graph::G
    pos::Int
    depth::Int
    pred::Dict{Int,Union{Int,Nothing}}
    level::Dict{Int,Int}
    low::Dict{Int,Int}
    bridges::Vector{Int}
end

FindBridgeState(graph, root) = FindBridgeState(
    graph, root, 0, Dict{Int,Union{Int,Nothing}}(root => nothing),
    Dict(root => 1), Dict(root => 1), Int[]
)


function yieldedge(state::FindBridgeState, edge)
    push!(state.bridges, edge)
end


function bridges(graph::UDGraph, root)
    state = FindBridgeState(graph, root)
    findbridge(state)
    return state.bridges
end


function findbridge(state::FindBridgeState)
    i = state.pos
    state.depth += 1
    state.level[i] = state.depth
    state.low[i] = state.depth
    for nbr in neighborkeys(state.graph, i)
        if !(nbr in keys(state.pred))
            # New node
            state.pred[nbr] = i
            state.pos = nbr
            findbridge(state)
            if state.low[nbr] == state.level[i]
                # Articulation point (bridgehead)
            elseif state.low[nbr] > state.level[i]
                # Articulation point
                yieldedge(state, neighbors(state.graph, i)[nbr])
            end
            state.low[i] = min(state.low[i], state.low[nbr])
        elseif nbr != state.pred[i]
            # Cycle found
            state.low[i] = min(state.low[i], state.level[nbr])
        end
    end
end
