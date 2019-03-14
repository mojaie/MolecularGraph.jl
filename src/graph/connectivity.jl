#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    connected_components,
    articulation_points, bridges, biconnected_components,
    two_edge_connected


struct ConnectedComponentState{G<:UDGraph}
    graph::G

    visited::Set{Int}
    remaining::Set{Int}

    components::Vector{Set{Int}}

    function ConnectedComponentState{G}(graph) where {G<:UDGraph}
        new(graph, Set(), nodekeys(graph), Set{Int}[])
    end
end


function dfs!(state::ConnectedComponentState, n)
    push!(state.visited, n)
    for nbr in neighborkeys(state.graph, n)
        if !(nbr in state.visited)
            dfs!(state, nbr)
        end
    end
end


function run!(state::ConnectedComponentState)
    while !isempty(state.remaining)
        dfs!(state, pop!(state.remaining))
        push!(state.components, copy(state.visited))
        setdiff!(state.remaining, state.visited)
        empty!(state.visited)
    end
end


"""
    connected_components(graph::UDGraph) -> Set{Set{Int}}

Compute connectivity and return sets of the connected components.
"""
function connected_components(graph::G) where {G<:UDGraph}
    state = ConnectedComponentState{G}(graph)
    run!(state)
    return state.components
end



struct BiconnectedState{G<:UDGraph}
    graph::G

    pred::Dict{Int,Int}
    level::Dict{Int,Int}
    low::Dict{Int,Int}
    compbuf::Set{Int}

    cutvertices::Set{Int}
    bridges::Set{Int}
    biconnected::Vector{Set{Int}}

    function BiconnectedState{G}(graph) where {G<:UDGraph}
        new(graph, Dict(), Dict(), Dict(), Set(), Set(), Set(), Set{Int}[])
    end
end


dfs!(state::BiconnectedState) = dfs!(state, 1, pop!(nodekeys(state.graph)))


function dfs!(state::BiconnectedState, depth::Int, n::Int)
    state.level[n] = depth
    state.low[n] = depth
    for (nbr, bond) in neighbors(state.graph, n)
        if n in keys(state.pred) && nbr == state.pred[n]
            continue # predecessor
        elseif !(nbr in keys(state.pred))
            # New node
            state.pred[nbr] = n
            push!(state.compbuf, n)
            dfs!(state, depth + 1, nbr)
            if state.low[nbr] >= state.level[n]
                # Articulation point
                if state.low[nbr] > state.level[n]
                    push!(state.bridges, bond) # except for bridgehead
                end
                push!(state.cutvertices, n)
                push!(state.biconnected, copy(state.compbuf))
                empty!(state.compbuf)
            end
            state.low[n] = min(state.low[n], state.low[nbr])
        else
            # Cycle found
            state.low[n] = min(state.low[n], state.level[nbr])
        end
    end
end


"""
    articulation_points(graph::UDGraph) -> Set{Int}

Compute biconnectivity and return articulation points.
"""
function articulation_points(graph::G) where {G<:UDGraph}
    state = BiconnectedState{G}(graph)
    dfs!(state)
    return state.cutvertices
end


"""
    bridges(graph::UDGraph) -> Set{Int}

Compute biconnectivity and return bridges.
"""
function bridges(graph::G) where {G<:UDGraph}
    state = BiconnectedState{G}(graph)
    dfs!(state)
    return state.bridges
end


"""
    bridges(graph::UDGraph) -> Set{Int}

Compute biconnectivity and return sets of the biconnected components.
"""
function biconnected_components(graph::G) where {G<:UDGraph}
    state = BiconnectedState{G}(graph)
    dfs!(state)
    return state.biconnected
end


"""
    two_edge_connected(graph::UDGraph) -> Set{Int}

Compute biconnectivity and return sets of the 2-edge connected components.
Isolated nodes will be filtered out.
"""
function two_edge_connected(graph::G) where {G<:UDGraph}
    brs = Int[]
    for conn in connected_components(graph)
        subg = nodesubgraph(graph, conn)
        state = BiconnectedState{G}(graph)
        dfs!(state)
        append!(brs, state.bridges)
    end
    bicomp = UDSubgraph(graph, nodekeys(graph), setdiff(edgekeys(graph), brs))
    return Set{Int}[c for c in connected_components(bicomp) if length(c) > 1]
end
