#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    @cache, clearcache!,
    AbstractGraph,
    UndirectedGraph, DirectedGraph, OrderedGraph, OrderedDiGraph,
    AbstractNode, AbstractEdge, UndirectedEdge, DirectedEdge,
    neighbors, outneighbors, inneighbors,
    findedgekey, getedge, hasedge, nodeattr, edgeattr,
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
abstract type OrderedGraph <: UndirectedGraph end
abstract type OrderedDiGraph <: DirectedGraph end

# TODO: use traits
# https://github.com/JuliaLang/julia/issues/2345



# Components

abstract type AbstractNode end
abstract type AbstractEdge end
abstract type UndirectedEdge <: AbstractEdge end
abstract type DirectedEdge <: AbstractEdge end



# Indexing

Base.getindex(graph::AbstractGraph, sym::Symbol) = eval(Expr(:call, sym, graph))
Base.getindex(
    graph::AbstractGraph, k1::Symbol, k2::Symbol, K::Symbol...
) = hcat(eval(Expr(:call, sym, graph)) for k in [k1, k2, K...])



# Cached properties

macro cache(ex)
    func = ex.args[1].args[1]
    dummy = gensym()
    ex.args[1].args[1] = dummy
    return quote
        $(esc(ex))
        Core.@__doc__ function $(esc(func))(graph)
            if !isdefined(graph, :cache)
                # Cache not available
                return $(esc(dummy))(graph)
            end
            symf = nameof($(esc(func)))
            if symf in keys(graph.cache)
                # Return cache
                return graph.cache[symf]
            end
            # Otherwise, set cache
            graph.cache[symf] = $(esc(dummy))(graph)
            return graph.cache[symf]
        end
    end
end


function clearcache!(graph)
    empty!(graph.cache)
    return
end



# Lookup

"""
    neighbors(graph, i) -> Dict{Int,Int}

Return the mapping of incident edges and adjacent nodes of node `i`.
If the graph is directed graph, both outneighbors and inneighbors are mapped.
"""
neighbors(graph::UndirectedGraph, i::Int) = graph.neighbormap[i]
neighbors(graph::DirectedGraph, i::Int
    ) = merge(outneighbors(graph, i), inneighbors(graph, i))


"""
    outneighbors(graph::DirectedGraph, i::Int) -> Dict{Int,Int}

Return the mapping of successor node keys and out edge keys connected to
the given node.
"""
outneighbors(graph::DirectedGraph, i::Int) = graph.outneighbormap[i]


"""
    inneighbors(graph::DirectedGraph, i::Int) -> Dict{Int,Int}

Return the mapping of predecessor node keys and in edge keys connected to
the given node.
"""
inneighbors(graph::DirectedGraph, i::Int) = graph.inneighbormap[i]


function findedgekey(graph::UndirectedGraph, u::Int, v::Int)
    for (inc, adj) in neighbors(graph, u)
        adj == v && return inc
    end
end

function findedgekey(graph::DirectedGraph, source::Int, target::Int)
    for (inc, adj) in outneighbors(graph, source)
        adj == target && return inc
    end
end


"""
    getedge(graph::AbstractGraph, i::Int) -> Tuple{Int,Int}

Return an edge tuple.
"""
getedge(graph::AbstractGraph, i::Int) = graph.edges[i]


"""
    hasedge(graph::UndirectedGraph, u::Int, v::Int) -> Bool
    hasedge(graph::DirectedGraph, source::Int, target::Int) -> Bool

Return whether the given two nodes are connected by at least one edge.
"""
hasedge(graph::AbstractGraph, u::Int, v::Int
    ) = findedgekey(graph, u, v) !== nothing


"""
    nodeattr(graph::AbstractGraph, i::Int) -> AbstractNode

Return the attribute object of node `i`.
"""
nodeattr(graph::AbstractGraph, i::Int) = graph.nodeattrs[i]


"""
    edgeattr(graph::AbstractGraph, i::Int) -> AbstractEdge

Return the attribute object of edge `i`.
"""
edgeattr(graph::AbstractGraph, i::Int) = graph.edgeattrs[i]

"""
    edgeattr(graph::AbstractGraph, u::Int, v::Int
        ) -> Union{AbstractEdge,Nothing}

Return the attribute object of an edge that connects `u` and `v`. If not found,
return nothing.
"""
function edgeattr(graph::AbstractGraph, u::Int, v::Int)
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
    nodeset(graph::Union{OrderedGraph,OrderedDiGraph}) -> Set{Int}

Return the set of node keys.
"""
nodeset(graph::Union{OrderedGraph,OrderedDiGraph}) = Set(1:nodecount(graph))


"""
    edgeset(graph::Union{OrderedGraph,OrderedDiGraph}) -> Set{Int}

Return the set of edge keys.
"""
edgeset(graph::Union{OrderedGraph,OrderedDiGraph}) = Set(1:edgecount(graph))


edgesiter(graph::Union{OrderedGraph,OrderedDiGraph}) = graph.edges
neighborsiter(graph::Union{OrderedGraph}) = graph.neighbormap


"""
    nodeattrs(graph::Union{OrderedGraph,OrderedDiGraph}) -> Vector{AbstractNode}

Return graph node attributes.
"""
nodeattrs(graph::Union{OrderedGraph,OrderedDiGraph}) = graph.nodeattrs


"""
    edgeattrs(graph::Union{OrderedGraph,OrderedDiGraph}) -> Vector{AbstractEdge}

Return graph edge attributes.
"""
edgeattrs(graph::Union{OrderedGraph,OrderedDiGraph}) = graph.edgeattrs



# Counting

"""
    nodecount(graph::AbstractGraph) -> Int

Return the number of graph nodes.
"""
nodecount(graph::UndirectedGraph) = length(graph.neighbormap)
nodecount(graph::DirectedGraph) = length(graph.outneighbormap)


"""
    edgecount(graph::AbstractGraph) -> Int

Return the number of graph edges.
"""
edgecount(graph::AbstractGraph) = length(graph.edges)


"""
    neighborcount(graph::AbstractGraph, n::Int) -> Int
    degree(graph::AbstractGraph, n::Int) -> Int

Return the number of adjacent nodes of the node 'n'.
"""
neighborcount(graph::AbstractGraph, n::Int) = length(neighbors(graph, n))
degree = neighborcount


"""
    outdegree(graph::DirectedGraph, n::Int) -> Int

Return the number of outneighbors of the node 'n'.
"""
outdegree(graph::DirectedGraph, n::Int) = length(outneighbors(graph, n))


"""
    indegree(graph::DirectedGraph, n::Int) -> Int

Return the number of inneighbors of the node 'n'.
"""
indegree(graph::DirectedGraph, n::Int) = length(inneighbors(graph, n))



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


"""
    addedge!(graph, u, v) -> Int
    addedge!(graph, u, v, attr) -> Int

Add new edge and return the edge index. If the edge attribute type is required,
specify the edge attribute object by `attr` keyword.
"""
function _addedge!(graph::UndirectedGraph, u::Int, v::Int)
    push!(graph.edges, (u, v))
    i = edgecount(graph)
    graph.neighbormap[u][i] = v
    graph.neighbormap[v][i] = u
    return i
end

function addedge!(graph::UndirectedGraph, u::Int, v::Int, attr::AbstractEdge)
    push!(graph.edgeattrs, attr)
    return _addedge!(graph, u, v)
end

function addedge!(graph::UndirectedGraph, u::Int, v::Int)
    isdefined(graph, :edgeattrs) && throw(ErrorException("edgeattr required"))
    return _addedge!(graph, u, v)
end

function _addedge!(graph::DirectedGraph, s::Int, t::Int)
    push!(graph.edges, (s, t))
    i = edgecount(graph)
    graph.outneighbormap[s][i] = t
    graph.inneighbormap[t][i] = s
    return i
end

function addedge!(graph::DirectedGraph, s::Int, t::Int, attr::AbstractEdge)
    push!(graph.edgeattrs, attr)
    return _addedge!(graph, s, t)
end

function addedge!(graph::DirectedGraph, s::Int, t::Int)
    isdefined(graph, :edgeattrs) && throw(ErrorException("edgeattr required"))
    return _addedge!(graph, s, t)
end


"""
    setnodeattr!(graph::AbstractGraph, i::Int, attr::AbstractNode)

Update the node attribute.
"""
function setnodeattr!(graph::AbstractGraph, i::Int, attr::AbstractNode)
    graph.nodeattrs[i] = attr
    return
end


"""
    setedgeattr!(graph::AbstractGraph, i::Int, attr::AbstractNode)

Update the edge attribute.
"""
function setedgeattr!(graph::AbstractGraph, i::Int, attr::AbstractEdge)
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
