#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    connectedcomponents, connectedmembership,
    cutvertices, bridges, biconnectedcomponents, biconnectedmembership,
    twoedgeconnectedcomponents, twoedgemembership


"""
    connectedcomponents(graph::UndirectedGraph) -> Vector{Set{Int}}

Compute connectivity and return sets of the connected components.
"""
function connectedcomponents(graph::UndirectedGraph)
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

connectedcomponents(view::SubgraphView) = connectedcomponents(view.graph)


"""
    connectedmembership(graph::OrderedGraph) -> Vector{Int}

Return connected component membership array.
"""
function connectedmembership(graph::OrderedGraph)
    mem = zeros(Int, nodecount(graph))
    for (i, conn) in enumerate(connectedcomponents(graph))
        for c in conn
            mem[c] = i
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
    return getproperty(state, sym)
end


"""
    cutvertices(graph::UndirectedGraph) -> Set{Int}

Compute biconnectivity and return cut vertices (articulation points).
"""
function cutvertices(graph::UndirectedGraph)
    return findbiconnected(graph, :cutvertices)
end

cutvertices(view::SubgraphView) = cutvertices(view.graph)


"""
    bridges(graph::UndirectedGraph) -> Set{Int}

Compute biconnectivity and return bridges.
"""
function bridges(graph::UndirectedGraph)
    return findbiconnected(graph, :bridges)
end

bridges(view::SubgraphView) = bridges(view.graph)


"""
    biconnected(graph::UndirectedGraph) -> Vector{Vector{Int}}

Compute sets of biconnected component edges.
"""
function biconnectedcomponents(graph::UndirectedGraph)
    return findbiconnected(graph, :biconnected)
end

biconnectedcomponents(view::SubgraphView) = biconnectedcomponents(view.graph)


function biconnectedmembership(graph::OrderedGraph)
    mem = zeros(Int, edgecount(graph))
    for (i, conn) in enumerate(biconnectedcomponents(graph))
        mem[conn] .= i
    end
    return mem
end

biconnectedmembership(view::SubgraphView) = biconnectedmembership(view.graph)


"""
    twoedgeconnectedcomponents(graph::UndirectedGraph) -> Vector{Set{Int}}

Compute sets of 2-edge connected component nodes.
"""
function twoedgeconnectedcomponents(graph::UndirectedGraph)
    cobr = setdiff(edgeset(graph), bridges(graph))
    subg = plaingraph(SubgraphView(graph, nodeset(graph), cobr))
    return connectedcomponents(subg)
end

twoedgeconnectedcomponents(view::SubgraphView
    ) = twoedgeconnectedcomponents(view.graph)


function twoedgemembership(graph::OrderedGraph)
    mem = zeros(Int, nodecount(graph))
    for (i, conn) in enumerate(twoedgeconnectedcomponents(graph))
        mem[conn] .= i
    end
    return mem
end

twoedgemembership(view::SubgraphView) = twoedgemembership(view.graph)
