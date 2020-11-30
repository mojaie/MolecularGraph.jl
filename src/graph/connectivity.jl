#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    connectedcomponents, connectedmembership,
    cutvertices, bridges, edgebiconnectedcomponents


"""
    connectedcomponents(graph::UndirectedGraph) -> Vector{Set{Int}}

Compute connectivity and return sets of the connected components.
"""
@cachefirst function connectedcomponents(graph::UndirectedGraph)
    nodes = nodeset(graph)
    components = Set{Int}[]
    while !isempty(nodes)
        root = pop!(nodes)
        tree = Set(keys(dfstree(adjacencies, graph, root)))
        push!(components, tree)
        setdiff!(nodes, tree)
    end
    return components
end

connectedcomponents(view::SubgraphView) = connectedcomponents(view.graph)


"""
    connectedmembership(graph::OrderedGraph) -> Vector{Int}

Return connected component membership array.
"""
@cachefirst function connectedmembership(graph::OrderedGraph)
    mem = zeros(Int, nodecount(graph))
    for (i, conn) in enumerate(connectedcomponents(graph))
        for n in conn
            mem[n] = i
        end
    end
    return mem
end

connectedmembership(view::SubgraphView) = connectedmembership(view.graph)


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
BiconnectedState(graph::UndirectedGraph) = BiconnectedState{typeof(graph)}(graph)


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


function findbiconnected(graph::UndirectedGraph, sym::Symbol)
    state = BiconnectedState(graph)
    nodes = nodeset(graph)
    while !isempty(nodes)
        dfs!(state, pop!(nodes))
        setdiff!(nodes, keys(state.level))
    end
    return getproperty(state, sym)
end


"""
    cutvertices(graph::UndirectedGraph) -> Set{Int}

Compute biconnectivity and return cut vertices (articulation points).
"""
@cachefirst function cutvertices(graph::UndirectedGraph)
    return findbiconnected(graph, :cutvertices)
end

cutvertices(view::SubgraphView) = cutvertices(view.graph)


"""
    bridges(graph::UndirectedGraph) -> Set{Int}

Compute biconnectivity and return bridges.
"""
@cachefirst function bridges(graph::UndirectedGraph)
    return findbiconnected(graph, :bridges)
end

bridges(view::SubgraphView) = bridges(view.graph)


"""
    edgebiconnectedcomponents(graph::UndirectedGraph) -> Vector{Vector{Int}}

Compute sets of biconnected component edges.
"""
@cachefirst function edgebiconnectedcomponents(graph::UndirectedGraph)
    return findbiconnected(graph, :biconnected)
end

edgebiconnectedcomponents(view::SubgraphView
    ) = edgebiconnectedcomponents(view.graph)
