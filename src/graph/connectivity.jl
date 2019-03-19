#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    connected_components, connected_membership,
    articulation_points, bridges,
    biconnected_components, biconnected_membership,
    two_edge_connected, two_edge_membership


struct ConnectedComponentState{G<:UDGraph}
    graph::G

    visited::Set{Int}
    remaining::Set{Int}

    components::Vector{Vector{Int}}

    function ConnectedComponentState{G}(graph) where {G<:UDGraph}
        new(graph, Set(), nodeset(graph), [])
    end
end


function dfs!(state::ConnectedComponentState, n)
    push!(state.visited, n)
    for nbr in neighborset(state.graph, n)
        if !(nbr in state.visited)
            dfs!(state, nbr)
        end
    end
end


function run!(state::ConnectedComponentState)
    while !isempty(state.remaining)
        dfs!(state, pop!(state.remaining))
        push!(state.components, collect(state.visited))
        setdiff!(state.remaining, state.visited)
        empty!(state.visited)
    end
end


"""
    connected_components(graph::UDGraph) -> Vector{Vector{Int}}

Compute connectivity and return sets of the connected components.
"""
function connected_components(graph::G) where {G<:UDGraph}
    hasproperty(graph, :connected) && return graph.property.connected
    state = ConnectedComponentState{G}(graph)
    run!(state)
    if isdefined(graph, :property)
        append!(graph.property.connected, state.components)
        return graph.property.connected
    end
    return state.components
end

function connected_membership(graph::UDGraph)
    hasproperty(graph, :connmembership) && return graph.property.connmembership
    mem = zeros(Int, nodecount(graph))
    for (i, conn) in enumerate(connected_components(graph))
        for c in conn
            mem[c] = i
        end
    end
    if isdefined(graph, :property)
        append!(graph.property.connmembership, mem)
        return graph.property.connmembership
    end
    return mem
end


struct BiconnectedState{G<:UDGraph}
    graph::G

    pred::Dict{Int,Int}
    level::Dict{Int,Int}
    low::Dict{Int,Int}
    compbuf::Vector{Int}

    cutvertices::Vector{Int}
    bridges::Vector{Int}
    biconnected::Vector{Vector{Int}}

    function BiconnectedState{G}(graph) where {G<:UDGraph}
        new(graph, Dict(), Dict(), Dict(), [], [], [], [])
    end
end


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


function biconnected(graph::G, sym::Symbol) where {G<:UDGraph}
    hasproperty(graph, sym) && return getproperty(graph.property, sym)
    state = BiconnectedState{G}(graph)
    nodes = nodeset(graph)
    while !isempty(nodes)
        dfs!(state, 1, pop!(nodes))
        setdiff!(nodes, keys(state.level))
    end
    if isdefined(graph, :property)
        append!(graph.property.biconnected, state.biconnected)
        append!(graph.property.cutvertices, state.cutvertices)
        append!(graph.property.bridges, state.bridges)
        return getproperty(graph.property, sym)
    else
        return getproperty(state, sym)
    end
end


"""
    articulation_points(graph::UDGraph) -> Set{Int}

Compute biconnectivity and return articulation points.
"""
articulation_points(graph::UDGraph) = biconnected(graph, :cutvertices)


"""
    bridges(graph::UDGraph) -> Set{Int}

Compute biconnectivity and return bridges.
"""
bridges(graph::UDGraph) = biconnected(graph, :bridges)


"""
    biconnected_components(graph::UDGraph) -> Vector{Set{Int}}

Compute biconnectivity and return sets of biconnected components.
"""
biconnected_components(graph::UDGraph) = biconnected(graph, :biconnected)


function biconnected_membership(graph::UDGraph)
    hasproperty(graph, :biconnmembership) && return graph.property.biconnmembership
    mem = zeros(Int, nodecount(graph))
    for (i, conn) in enumerate(biconnected_components(graph))
        for c in conn
            mem[c] = i
        end
    end
    if isdefined(graph, :property)
        append!(graph.property.biconnmembership, mem)
        return graph.property.biconnmembership
    end
    return mem
end


"""
    two_edge_connected(graph::UDGraph) -> Set{Int}

Compute biconnectivity and return sets of the 2-edge connected components.
Isolated nodes will be filtered out.
"""
function two_edge_connected(graph::UDGraph)
    hasproperty(graph, :twoedge) && return graph.property.twoedge
    cobr = setdiff(edgeset(graph), bridges(graph))
    comp = connected_components(edgesubgraph(graph, cobr))
    if isdefined(graph, :property)
        append!(graph.property.twoedge, comp)
        return graph.property.twoedge
    end
    return comp
end


function two_edge_membership(graph::UDGraph)
    hasproperty(graph, :twoedgemembership) && return graph.property.twoedgemembership
    mem = zeros(Int, nodecount(graph))
    for (i, conn) in enumerate(two_edge_connected(graph))
        for c in conn
            mem[c] = i
        end
    end
    if isdefined(graph, :property)
        append!(graph.property.twoedgemembership, mem)
        return graph.property.twoedgemembership
    end
    return mem
end
