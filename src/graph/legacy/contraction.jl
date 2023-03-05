#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    ContrationView, edgecontraction


struct ContractionView{T<:UndirectedGraph} <: UndirectedGraph
    graph::T
    collapsed::Set{Int}
    vnode::Int
    vneighbormap::Dict{Int, Dict{Int,Int}}
end


"""
    edgecontraction(graph::UndirectedGraph, nodes) -> ContractionView

Generate edge-contraction view.
"""
function edgecontraction(graph::UndirectedGraph, nodes)
    # Determine outgoing edges to be contracted
    incs = Set{Int}()
    for n in nodes
        union!(incs, incidences(graph, n))
    end
    inedges = edgeset(nodesubgraph(graph, nodes))
    outedges = setdiff(incs, inedges)

    # Virtual node
    vnode = maximum(nodeset(graph)) + 1

    # Neighbormap
    vnbrmap = Dict(vnode => Dict{Int,Int}())
    for e in outedges
        (u, v) = getedge(graph, e)
        if u in nodes
            if !(v in keys(vnbrmap))
                vnbrmap[v] = copy(neighbors(graph, v))
            end
            vnbrmap[v][e] = vnode
            vnbrmap[vnode][e] = v
        elseif v in nodes
            if !(u in keys(vnbrmap))
                vnbrmap[u] = copy(neighbors(graph, u))
            end
            vnbrmap[u][e] = vnode
            vnbrmap[vnode][e] = u
        else
            @assert false
        end
    end
    return ContractionView(graph, Set(nodes), vnode, vnbrmap)
end


function neighbors(view::ContractionView, idx::Int)
    idx in view.collapsed && throw(KeyError(idx))
    if idx in keys(view.vneighbormap)
        return view.vneighbormap[idx]
    else
        return neighbors(view.graph, idx)
    end
end


function getedge(view::ContractionView, idx::Int)
    (u, v) = getedge(view.graph, idx)
    if u in view.collapsed
        return (view.vnode, v)
    elseif v in view.collapsed
        return (u, view.vnode)
    else
        return (u, v)
    end
end


function hasedge(view::ContractionView, u::Int, v::Int)
    u in view.collapsed && throw(KeyError(idx))
    v in view.collapsed && throw(KeyError(idx))
    return findedgekey(view.graph, u, v) !== nothing
end


nodeset(view::ContractionView) = union(setdiff(nodeset(view.graph), view.collapsed), view.vnode)
edgeset(view::ContractionView) = edgeset(view.graph)

nodecount(view::ContractionView) = nodecount(view.graph) - length(view.collapsed) + 1
edgecount(view::ContractionView) = edgecount(view.graph)
