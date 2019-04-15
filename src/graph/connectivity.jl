#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    connected_components, connected_membership,
    cutvertices, bridges, biconnected_components, biconnected_membership,
    two_edge_connected, two_edge_membership


"""
    connected_components(graph::UndirectedGraph) -> Vector{Set{Int}}

Compute connectivity and return sets of the connected components.
"""
@cache function connected_components(graph::UndirectedGraph)
    nodes = nodeset(graph)
    components = Set{Int}[]
    while !isempty(nodes)
        root = pop!(nodes)
        tree = dfstree(adjacencies, graph, root)
        push!(components, union(keys(tree), [root]))
        setdiff!(nodes, keys(tree))
    end
    return components
end

connected_components(view::SubgraphView) = connected_components(view.graph)


"""
    connected_membership(graph::OrderedGraph) -> Vector{Int}

Return connected component membership array.
"""
@cache function connected_membership(graph::OrderedGraph)
    mem = zeros(Int, nodecount(graph))
    for (i, conn) in enumerate(connected_components(graph))
        for c in conn
            mem[c] = i
        end
    end
    return mem
end

connected_membership(view::SubgraphView) = connected_membership(view.graph)


struct BiconnectedState{T<:UndirectedGraph}
    graph::T

    pred::Dict{Int,Int}
    level::Dict{Int,Int}
    low::Dict{Int,Int}

    cutvertices::Set{Int}
    bridges::Vector{Int}
    biconnected::Vector{Set{Int}} # biconnected edges

    function BiconnectedState{T}(graph) where {T<:UndirectedGraph}
        new(graph, Dict(), Dict(), Dict(), Set(), [], [])
    end
end


function dfs!(state::BiconnectedState, n::Int)
    state.pred[n] = n
    dfs!(state, 1, n)
end

function dfs!(state::BiconnectedState, depth::Int, n::Int)
    state.level[n] = depth
    state.low[n] = depth
    childcnt = 0
    compbuf = Set{Int}()
    for (ninc, nadj) in neighbors(state.graph, n)
        if !haskey(state.level, nadj)
            # New node
            childcnt += 1
            state.pred[nadj] = n
            comp = dfs!(state, depth + 1, nadj)
            push!(comp, ninc)
            if state.low[nadj] >= state.level[n]
                # Articulation point
                if state.low[nadj] > state.level[n]
                    push!(state.bridges, ninc) # except for bridgehead
                end
                push!(state.biconnected, comp)
                push!(state.cutvertices, n)
            else
                union!(compbuf, comp)
            end
            state.low[n] = min(state.low[n], state.low[nadj])
        elseif state.pred[n] != nadj
            # Cycle found
            state.low[n] = min(state.low[n], state.level[nadj])
            state.pred[nadj] = n
            push!(compbuf, ninc)
        end
    end
    if depth == 1 && childcnt < 2
        delete!(state.cutvertices, n)
    end
    return compbuf
end


function findbiconnected(graph::T, sym::Symbol) where {T<:UndirectedGraph}
    state = BiconnectedState{T}(graph)
    nodes = nodeset(graph)
    while !isempty(nodes)
        dfs!(state, pop!(nodes))
        setdiff!(nodes, keys(state.level))
    end
    if isdefined(graph, :cache)
        graph.cache[:biconnected_components] = state.biconnected
        graph.cache[:cutvertices] = state.cutvertices
        graph.cache[:bridges] = state.bridges
        return graph.cache[sym]
    else
        return getproperty(state, sym)
    end
end


"""
    cutvertices(graph::UndirectedGraph) -> Set{Int}

Compute biconnectivity and return cut vertices (articulation points).
"""
@cache function cutvertices(graph::UndirectedGraph)
    return findbiconnected(graph, :cutvertices)
end

cutvertices(view::SubgraphView) = cutvertices(view.graph)


"""
    bridges(graph::UndirectedGraph) -> Set{Int}

Compute biconnectivity and return bridges.
"""
@cache function bridges(graph::UndirectedGraph)
    return findbiconnected(graph, :bridges)
end

bridges(view::SubgraphView) = bridges(view.graph)


"""
    biconnected(graph::UndirectedGraph) -> Vector{Vector{Int}}

Compute sets of biconnected component edges.
"""
@cache function biconnected_components(graph::UndirectedGraph)
    return findbiconnected(graph, :biconnected_components)
end

biconnected_components(view::SubgraphView) = biconnected_components(view.graph)


@cache function biconnected_membership(graph::OrderedGraph)
    mem = zeros(Int, edgecount(graph))
    for (i, conn) in enumerate(biconnected_components(graph))
        mem[conn] .= i
    end
    return mem
end

biconnected_membership(view::SubgraphView) = biconnected_membership(view.graph)


"""
    two_edge_connected(graph::UndirectedGraph) -> Vector{Set{Int}}

Compute sets of 2-edge connected component nodes.
"""
@cache function two_edge_connected(graph::UndirectedGraph)
    cobr = setdiff(edgeset(graph), bridges(graph))
    subg = plaingraph(SubgraphView(graph, nodeset(graph), cobr))
    return connected_components(subg)
end

two_edge_connected(view::SubgraphView) = two_edge_connected(view.graph)


@cache function two_edge_membership(graph::OrderedGraph)
    mem = zeros(Int, nodecount(graph))
    for (i, conn) in enumerate(two_edge_connected(graph))
        mem[conn] .= i
    end
    return mem
end

two_edge_membership(view::SubgraphView) = two_edge_membership(view.graph)
