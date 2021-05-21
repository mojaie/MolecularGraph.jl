#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    AbstractGraph,
    UndirectedGraph, DirectedGraph, HyperGraph,
    OrderedGraph, OrderedDiGraph, OrderedHyperGraph,
    AbstractNode, AbstractEdge, UndirectedEdge, DirectedEdge,
    @cachefirst, @cache, hascache, clearcache!, clone,
    neighbors, outneighbors, inneighbors,
    findedgekey, findalledgekeys, getedge, hasedge, nodeattr, edgeattr,
    adjacencies, successors, predecessors,
    incidences, outincidences, in_incidences,
    nodeset, edgeset, edgesiter, neighborsiter, nodeattrs, edgeattrs,
    nodecount, edgecount, neighborcount, degree, indegree, outdegree,
    addnode!, addedge!, setnodeattr!, setedgeattr!,
    unlinknodes, unlinkedges,
    nodeattrtype, edgeattrtype



# Graph types

abstract type AbstractGraph end
abstract type UndirectedGraph <: AbstractGraph end
abstract type DirectedGraph <: AbstractGraph end
abstract type HyperGraph <: AbstractGraph end
abstract type OrderedGraph <: UndirectedGraph end
abstract type OrderedDiGraph <: DirectedGraph end
abstract type OrderedHyperGraph <: HyperGraph end


# Node and edge attributes

abstract type AbstractNode end
abstract type AbstractEdge end
abstract type UndirectedEdge <: AbstractEdge end # TODO: unnecessary?
abstract type DirectedEdge <: AbstractEdge end # TODO: unnecessary?



"""
    @cachefirst expression

A macro for molecular descriptor array cache mechanism.
"""
macro cachefirst(ex)
    func = ex.args[1].args[1]
    dummy = gensym()
    ex.args[1].args[1] = dummy
    return esc(quote
        $ex
        Core.@__doc__ function $func(graph; kwargs...)
            if isdefined(graph, :cache) && isempty(kwargs)
                # Return cache
                symf = nameof($func)
                symf in keys(graph.cache) && return graph.cache[symf]
            end
            return $dummy(graph; kwargs...)
        end
    end)
end


"""
    @cache expression

A macro that enables to force recalculate and set cache of descriptor arrays.
"""
macro cache(ex)
    func = ex.args[1]
    graph = ex.args[2]
    return esc(quote
        $graph.cache[nameof($func)] = $func($graph)
    end)
end


hascache(graph, key) = isdefined(graph, :cache) && haskey(graph.cache, key)


"""
    clearcache!(graph::AbstractGraph) -> nothing

Clear calculated property caches.

Calling `clearcache!` is recommended when the graph nodes/edges are added, removed or reindexed. You can [`Graph.clone`](@ref) the graph before destructive operation instead.
"""
function clearcache!(graph::AbstractGraph)
    empty!(graph.cache)
    return
end


"""
    clone(graph::AbstractGraph) -> AbstractGraph

Return deep copy of the graph.
"""
function clone end



# Lookup

"""
    neighbors(graph, i) -> Dict

Return the mapping of incident edges and adjacent nodes of node `i`.
If the graph is directed graph, both outneighbors and inneighbors are mapped.
"""
neighbors(graph::UndirectedGraph, i) = graph.neighbormap[i]
neighbors(graph::DirectedGraph, i) = merge(outneighbors(graph, i), inneighbors(graph, i))

neighbors(graph::HyperGraph, i) = Dict(inc => graph.edges[inc] for inc in graph.incidences[i])



"""
    outneighbors(graph::DirectedGraph, i) -> Dict

Return the mapping of successor node keys and out edge keys connected to
the given node.
"""
outneighbors(graph::DirectedGraph, i) = graph.outneighbormap[i]


"""
    inneighbors(graph::DirectedGraph, i) -> Dict

Return the mapping of predecessor node keys and in edge keys connected to
the given node.
"""
inneighbors(graph::DirectedGraph, i) = graph.inneighbormap[i]


function findedgekey(graph::UndirectedGraph, u, v)
    for (inc, adj) in neighbors(graph, u)
        adj == v && return inc
    end
end

function findedgekey(graph::DirectedGraph, source, target)
    for (inc, adj) in outneighbors(graph, source)
        adj == target && return inc
    end
end

function findedgekey(graph::HyperGraph, u, v)
    for (inc, adj) in neighbors(graph, u)
        v in adj && return inc
    end
end


function findalledgekeys(graph::UndirectedGraph, u, v)
    hasmultiedge(graph) || return [findedgekey(graph, u, v)]
    keys = Int[]
    for (inc, adj) in neighbors(graph, u)
        adj == v && push!(keys, inc)
    end
    return keys
end

function findalledgekeys(graph::DirectedGraph, source, target)
    hasmultiedge(graph) || return [findedgekey(graph, source, target)]
    keys = Int[]
    for (inc, adj) in outneighbors(graph, source)
        adj == target && push!(keys, inc)
    end
    return keys
end

function findalledgekeys(graph::HyperGraph, u, v)
    hasmultiedge(graph) || return [findedgekey(graph, u, v)]
    keys = Int[]
    for (inc, adj) in neighbors(graph, u)
        v in adj && push!(keys, inc)
    end
    return keys
end



"""
    getedge(graph::AbstractGraph, i) -> Any

Return an edge.
"""
getedge(graph::AbstractGraph, i) = graph.edges[i]



"""
    hasedge(graph::AbstractGraph, u, v) -> Bool

Return whether the given two nodes are connected by at least one edge.
"""
hasedge(graph::AbstractGraph, u, v) = findedgekey(graph, u, v) !== nothing


"""
    nodeattr(graph::AbstractGraph, i) -> AbstractNode

Return the attribute object of node `i`.
"""
nodeattr(graph::AbstractGraph, i) = graph.nodeattrs[i]


"""
    edgeattr(graph::AbstractGraph, i) -> AbstractEdge

Return the attribute object of edge `i`.
"""
edgeattr(graph::AbstractGraph, i) = graph.edgeattrs[i]

"""
    edgeattr(graph::AbstractGraph, u, v) -> Union{AbstractEdge,Nothing}

Return the attribute object of an edge that connects `u` and `v`. If not found,
return nothing.
"""
function edgeattr(graph::AbstractGraph, u, v)
    k = findedgekey(graph, u, v)
    return k === nothing ? nothing : graph.edgeattrs[k]
end



# Collection

adjacencies(g::AbstractGraph, i) = Set(values(neighbors(g, i)))
successors(g::DirectedGraph, i) = Set(values(outneighbors(g, i)))
predecessors(g::DirectedGraph, i) = Set(values(inneighbors(g, i)))

incidences(g::AbstractGraph, i) = Set(keys(neighbors(g, i)))
outincidences(g::DirectedGraph, i) = Set(keys(outneighbors(g, i)))
in_incidences(g::DirectedGraph, i) = Set(keys(inneighbors(g, i)))



"""
    nodeset(graph::Union{OrderedGraph,OrderedDiGraph,OrderedHyperGraph}) -> Set{Int}

Return the set of node keys.
"""
nodeset(graph::Union{OrderedGraph,OrderedDiGraph,OrderedHyperGraph}) = Set(1:nodecount(graph))


"""
    edgeset(graph::Union{OrderedGraph,OrderedDiGraph,OrderedHyperGraph}) -> Set{Int}

Return the set of edge keys.
"""
edgeset(graph::Union{OrderedGraph,OrderedDiGraph,OrderedHyperGraph}) = Set(1:edgecount(graph))


edgesiter(graph::Union{OrderedGraph,OrderedDiGraph,OrderedHyperGraph}) = graph.edges

neighborsiter(graph::OrderedGraph) = graph.neighbormap


"""
    nodeattrs(graph::Union{OrderedGraph,OrderedDiGraph,OrderedHyperGraph}) -> Vector{AbstractNode}

Return graph node attributes.
"""
nodeattrs(graph::Union{OrderedGraph,OrderedDiGraph,OrderedHyperGraph}) = graph.nodeattrs


"""
    edgeattrs(graph::Union{OrderedGraph,OrderedDiGraph,OrderedHyperGraph}) -> Vector{AbstractEdge}

Return graph edge attributes.
"""
edgeattrs(graph::Union{OrderedGraph,OrderedDiGraph,OrderedHyperGraph}) = graph.edgeattrs



# Counting

"""
    nodecount(graph::AbstractGraph) -> Int

Return the number of graph nodes.
"""
nodecount(graph::UndirectedGraph) = length(graph.neighbormap)
nodecount(graph::DirectedGraph) = length(graph.outneighbormap)
nodecount(graph::HyperGraph) = length(graph.incidences)


"""
    edgecount(graph::AbstractGraph) -> Int

Return the number of graph edges.
"""
edgecount(graph::AbstractGraph) = length(graph.edges)


"""
    neighborcount(graph::AbstractGraph, n) -> Int
    degree(graph::AbstractGraph, n) -> Int

Return the number of adjacent nodes of the node 'n'.
"""
neighborcount(graph::AbstractGraph, n) = length(neighbors(graph, n))
degree = neighborcount


"""
    outdegree(graph::DirectedGraph, n) -> Int

Return the number of outneighbors of the node 'n'.
"""
outdegree(graph::DirectedGraph, n) = length(outneighbors(graph, n))


"""
    indegree(graph::DirectedGraph, n) -> Int

Return the number of inneighbors of the node 'n'.
"""
indegree(graph::DirectedGraph, n) = length(inneighbors(graph, n))



# Editing

"""
    addnode!(graph) -> Int
    addnode!(graph, attr) -> Int

Add new node and return the node index. If the node attribute type is required,
specify the node attribute object by `node` keyword.
"""
function addnode!(graph::UndirectedGraph, attr::AbstractNode)
    push!(graph.neighbormap, Dict())
    push!(graph.nodeattrs, attr)
    return length(graph.neighbormap)
end

function addnode!(graph::UndirectedGraph)
    isdefined(graph, :nodeattrs) && throw(ErrorException("nodeattr required"))
    push!(graph.neighbormap, Dict())
    return length(graph.neighbormap)
end

function addnode!(graph::DirectedGraph, attr::AbstractNode)
    push!(graph.outneighbormap, Dict())
    push!(graph.inneighbormap, Dict())
    push!(graph.nodeattrs, attr)
    return length(graph.outneighbormap)
end

function addnode!(graph::DirectedGraph)
    isdefined(graph, :nodeattrs) && throw(ErrorException("nodeattr required"))
    push!(graph.outneighbormap, Dict())
    push!(graph.inneighbormap, Dict())
    return length(graph.outneighbormap)
end

function addnode!(graph::HyperGraph, attr::AbstractNode)
    push!(graph.incidences, Set())
    push!(graph.nodeattrs, attr)
    return length(graph.incidences)
end

function addnode!(graph::HyperGraph)
    isdefined(graph, :nodeattrs) && throw(ErrorException("nodeattr required"))
    push!(graph.incidences, Set())
    return length(graph.incidences)
end



"""
    addedge!(graph, u, v) -> Int
    addedge!(graph, u, v, attr) -> Int

Add new edge and return the edge index. If the edge attribute type is required,
specify the edge attribute object by `attr` keyword.
"""
function _addedge!(graph::UndirectedGraph, u, v)
    push!(graph.edges, (u, v))
    i = edgecount(graph)
    graph.neighbormap[u][i] = v
    graph.neighbormap[v][i] = u
    return i
end

function addedge!(graph::UndirectedGraph, u, v, attr::AbstractEdge)
    push!(graph.edgeattrs, attr)
    return _addedge!(graph, u, v)
end

function addedge!(graph::UndirectedGraph, u, v)
    isdefined(graph, :edgeattrs) && throw(ErrorException("edgeattr required"))
    return _addedge!(graph, u, v)
end

function _addedge!(graph::DirectedGraph, s, t)
    push!(graph.edges, (s, t))
    i = edgecount(graph)
    graph.outneighbormap[s][i] = t
    graph.inneighbormap[t][i] = s
    return i
end

function addedge!(graph::DirectedGraph, s, t, attr::AbstractEdge)
    push!(graph.edgeattrs, attr)
    return _addedge!(graph, s, t)
end

function addedge!(graph::DirectedGraph, s, t)
    isdefined(graph, :edgeattrs) && throw(ErrorException("edgeattr required"))
    return _addedge!(graph, s, t)
end


function _addedge!(graph::HyperGraph, edge)
    push!(graph.edges, edge)
    i = edgecount(graph)
    for n in edge
        push!(graph.incidences[n], i)
    end
    return i
end

function addedge!(graph::HyperGraph, edge, attr::AbstractEdge)
    push!(graph.edgeattrs, attr)
    return _addedge!(graph, edge)
end

function addedge!(graph::HyperGraph, edge)
    isdefined(graph, :edgeattrs) && throw(ErrorException("edgeattr required"))
    return _addedge!(graph, edge)
end



"""
    setnodeattr!(graph::AbstractGraph, i, attr::AbstractNode)

Update the node attribute.
"""
function setnodeattr!(graph::AbstractGraph, i, attr::AbstractNode)
    graph.nodeattrs[i] = attr
    return
end


"""
    setedgeattr!(graph::AbstractGraph, i, attr::AbstractNode)

Update the edge attribute.
"""
function setedgeattr!(graph::AbstractGraph, i, attr::AbstractEdge)
    graph.edgeattrs[i] = attr
    return
end


"""
    unlinknodes(graph, nodes)

Delete given nodes and its incident edges from the graph.
"""
function unlinknodes end


"""
    unlinkedges(graph, edges)

Delete given edges from the graph.
"""
function unlinkedges end



# Types

"""
    nodeattrtype(graph) -> Type

Return node attribute type.
"""
nodeattrtype(mol::AbstractGraph) = eltype(mol.nodeattrs)

"""
    edgeattrtype(graph) -> Type

Return edge attribute type.
"""
edgeattrtype(mol::AbstractGraph) = eltype(mol.edgeattrs)
