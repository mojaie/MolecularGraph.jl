#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

import Graphs:
    AbstractGraph, edgetype, nv, vertices, ne, edges,
    has_vertex, has_edge, outneighbors

export
    AbstractMolGraph, OrderedMolGraph, AbstractReaction,
    to_dict, edge_rank


abstract type AbstractMolGraph{T} <: AbstractGraph{T} end
abstract type OrderedMolGraph{T} <: AbstractMolGraph{T} end

abstract type AbstractReaction{T<:AbstractMolGraph} end


Base.copy(g::AbstractMolGraph) = deepcopy(mol)
Base.eltype(g::AbstractMolGraph) = eltype(typeof(g))
edgetype(g::AbstractMolGraph) = edgetype(typeof(g))

nv(g::AbstractMolGraph) = nv(g.graph)
vertices(g::AbstractMolGraph) = vertices(g.graph)
ne(g::AbstractMolGraph) = ne(g.graph)
edges(g::AbstractMolGraph) = edges(g.graph)

has_vertex(g::AbstractMolGraph, x::Integer) = has_vertex(g.graph, x)
has_edge(g::AbstractMolGraph, s::Integer, d::Integer) = has_edge(g.graph, s, d)

outneighbors(g::AbstractMolGraph, v::Integer) = outneighbors(g.graph, v)


function edge_rank(g::SimpleGraph, u::Integer, v::Integer)
    u, v = u < v ? (u, v) : (v, u)
    i = zero(u)
    cnt = 0
    @inbounds while i < u
        i += one(u)
        for j in g.fadjlist[i]
            if j > i
                cnt += 1
                j == v && break
            end
        end
    end
    return cnt
end

edge_rank(g::SimpleGraph, e::Edge) = edge_rank(g, src(e), dst(e))
