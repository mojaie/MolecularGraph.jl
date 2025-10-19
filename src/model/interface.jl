#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

"""
    AbstractMolGraph{T<:Integer} <: Graphs.AbstractGraph{T}

The base class of molecular graphs
"""
abstract type AbstractMolGraph{T<:Integer} <: Graphs.AbstractGraph{T} end


# Graphs.jl common interface

Base.eltype(g::AbstractMolGraph) = eltype(g.graph)
# https://github.com/JuliaGraphs/Graphs.jl/pull/144
Base.eltype(::Type{<:AbstractMolGraph{T}}) where T = T
Base.zero(::Type{G}) where G <: AbstractMolGraph = G()
Graphs.edgetype(g::AbstractMolGraph) = edgetype(g.graph)
Graphs.edgetype(::Type{<:AbstractMolGraph{T}}) where T = Edge{T}
Graphs.nv(g::AbstractMolGraph) = nv(g.graph)
Graphs.ne(g::AbstractMolGraph) = ne(g.graph)
Graphs.vertices(g::AbstractMolGraph) = vertices(g.graph)
Graphs.edges(g::AbstractMolGraph) = edges(g.graph)
Graphs.is_directed(::Type{<:AbstractMolGraph}) = false
Graphs.has_vertex(g::AbstractMolGraph, x::Integer) = has_vertex(g.graph, x)
Graphs.has_edge(g::AbstractMolGraph, s::Integer, d::Integer) = has_edge(g.graph, s, d)
Graphs.inneighbors(g::AbstractMolGraph, v::Integer) = inneighbors(g.graph, v)
Graphs.outneighbors(g::AbstractMolGraph, v::Integer) = outneighbors(g.graph, v)


"""
    vproptype(::Type{AbstractMolGraph}) -> Type
    vproptype(mol::AbstractMolGraph) -> Type

Return the type of vertex properties
"""
vproptype(::Type{T}) where T<:AbstractMolGraph = error("$(T) has no vertex property field")
vproptype(mol::T) where T<:AbstractMolGraph = vproptype(T)


"""
    eproptype(::Type{AbstractMolGraph}) -> Type
    eproptype(mol::AbstractMolGraph) -> Type

Return the type of edge properties
"""
eproptype(::Type{T}) where T<:AbstractMolGraph = error("$(T) has no edge property field")
eproptype(mol::T) where T<:AbstractMolGraph = eproptype(T)


edge_rank(mol::AbstractMolGraph) = edge_rank(mol.graph)



"""
    AbstractReaction

The base class of reactions
"""
abstract type AbstractReaction end



## Serializable integer keys

struct VertexKey{T<:Integer}
    key::T
end

VertexKey(x::VertexKey) = x
VertexKey{T}(x::VertexKey) where T <: Integer = VertexKey{T}(x.key)

Base.:(==)(x::VertexKey, y::VertexKey) = x.key == y.key
Base.:(==)(x::VertexKey{T}, y::T) where T = x.key == y  # required for getindex
Base.:(==)(x::T, y::VertexKey{T}) where T = x == y.key  # required for getindex
Base.:(+)(x::VertexKey{T}, y::T) where T = VertexKey{T}(x.key + y)
Base.hash(x::VertexKey, h::UInt) = hash(x.key, h)
Base.convert(::Type{VertexKey{T}}, x::T) where T = VertexKey{T}(x)
StructUtils.lowerkey(::JSON.JSONStyle, x::VertexKey{T}) where T = string(x.key)
StructUtils.liftkey(::Type{VertexKey{T}}, x::String) where T = VertexKey{T}(parse(T, x))


struct EdgeKey{T<:Integer}
    key::Edge{T}
end

EdgeKey(x::EdgeKey) = x
EdgeKey{T}(x::EdgeKey) where T <: Integer = EdgeKey{T}(x.key)

Base.:(==)(x::EdgeKey, y::EdgeKey) = x.key == y.key
Base.:(==)(x::EdgeKey{T}, y::Edge{T}) where T = x.key == y  # required for getindex
Base.:(==)(x::Edge{T}, y::EdgeKey{T}) where T = x == y.key  # required for getindex
Base.hash(x::EdgeKey, h::UInt) = hash(x.key, h)
Base.convert(::Type{EdgeKey{T}}, x::Edge{T}) where T = EdgeKey{T}(x)
StructUtils.lowerkey(::JSON.JSONStyle, x::EdgeKey{T}) where T = "$(x.key.src)_$(x.key.dst)"
StructUtils.liftkey(::Type{EdgeKey{T}}, x::String) where T = EdgeKey{T}(Edge([parse(T, s) for s in split(x, '_')]...))


"""
    AbstractProperty

The base class of graph-level property of the molecular model.
"""
abstract type AbstractProperty end

function Base.:(==)(g::T, h::T) where T <: AbstractProperty
    for sym in fieldnames(typeof(g))
        getproperty(g, sym) == getproperty(h, sym) || return false
    end
    return true
end


"""
    SimpleMolProperty{T<:Integer} <: AbstractProperty

The base class of graph-level property for SimpleMolGraph.
"""
abstract type SimpleMolProperty{T<:Integer} <: AbstractProperty end

Base.eltype(::Type{<:SimpleMolProperty{T}}) where T = T


"""
    AbstractState

The base class of molecular model states.
"""
abstract type AbstractState end


"""
    AbstractQueryNode

The base class of query node.
"""
abstract type AbstractQueryNode end

Base.getindex(elem::AbstractQueryNode, prop::Symbol) = getproperty(elem, prop)

Base.:(==)(x::T, y::T) where T <: AbstractQueryNode = all(
    [getfield(x, f) == getfield(y, f) for f in fieldnames(T)])

function Base.hash(elem::T, h::UInt) where T <: AbstractQueryNode
    for name in fieldnames(T)
        val = getfield(elem, name)
        h = hash(val, h)
    end
    return h
end


"""
    AbstractElement

The base class of vertex properties (atom) and edge properties (bond).
"""
abstract type AbstractElement end

Base.getindex(elem::AbstractElement, prop::Symbol) = getproperty(elem, prop)

Base.:(==)(x::T, y::T) where T <: AbstractElement = all(
    [getfield(x, f) == getfield(y, f) for f in fieldnames(T)])

function Base.hash(elem::T, h::UInt) where T <: AbstractElement
    for name in fieldnames(T)
        val = getfield(elem, name)
        h = hash(val, h)
    end
    return h
end



"""
    AbstractAtom <: AbstractElement

The base class of vertex properties (atom).
"""
abstract type AbstractAtom <: AbstractElement end

has_isaromatic(::Type{T}) where T <: AbstractAtom = false
has_mol(::Type{T}) where T <: AbstractAtom = false
has_formula(::Type{T}) where T <: AbstractAtom = false
has_hydrogens(::Type{T}) where T <: AbstractAtom = false
has_label(::Type{T}) where T <: AbstractAtom = false

"""
    atom_number(atom::AbstractAtom) -> Int

Return an atomic number of the given atom or the atomic symbol.
"""
atom_number(atom::AbstractAtom) = error("atom_number is not implemented for this atom type")


"""
    atom_symbol(atom::AbstractAtom) -> Symbol

Return an atomic symbol of the given atom or the atomic number.
"""
atom_symbol(atom::AbstractAtom) = error("atom_symbol is not implemented for this atom type")


"""
    atom_charge(atom::AbstractAtom) -> Int

Return atomic charge of the given atom.
"""
atom_charge(atom::AbstractAtom) = error("atom_charge is not implemented for this atom type")


"""
    multiplicity(atom::AbstractAtom) -> Int

Return multiplicity (num of radicals + 1) of the given atom.

This is experimental feature - free radical chemistry is still not introduced to this library.
This does nothing for now, but for example, you can set multiplicity=2 to molecular oxygens manually.
"""
multiplicity(atom::AbstractAtom) = error("multiplicity is not implemented for this atom type")



abstract type StandardAtom <: AbstractAtom end



"""
    AbstractBond <: AbstractElement

The base class of edge properties (bond).
"""
abstract type AbstractBond <: AbstractElement end

has_submap(::Type{T}) where T <: AbstractBond = false


"""
    bond_order(bond::AbstractBond) -> Int

Return bond order of the given bond.
"""
bond_order(bond::AbstractBond) = error("bond_order is not implemented for this bond type")


abstract type StandardBond <: AbstractBond end



"""
    QueryTree{T<:Integer,U<:AbstractQueryNode} <: AbstractElement

The base class of molecular query trees.
"""
abstract type QueryTree{T<:Integer,U<:AbstractQueryNode} <: AbstractElement end

Base.eltype(::Type{<:QueryTree{T,U}}) where {T,U} = T
Base.eltype(qtree::T) where T<:QueryTree = eltype(T)
vproptype(::Type{<:QueryTree{T,U}}) where {T,U} = U
vproptype(qtree::T) where T<:QueryTree = vproptype(T)



"""
    Base.getindex(mol::AbstractMolGraph, v::Integer) -> AbstractElement
    Base.getindex(mol::AbstractMolGraph, e::Edge) -> AbstractElement
    Base.getindex(mol::AbstractMolGraph, u::Integer, v::Integer) -> AbstractElement

Return properties (vertex or edge attributes).
"""
Base.getindex(mol::AbstractMolGraph, v::Integer) = mol.vprops[v]
Base.getindex(mol::AbstractMolGraph, e::Edge) = mol.eprops[e]
Base.getindex(mol::AbstractMolGraph, gprop::Symbol) = getproperty(mol.gprops, gprop)

"""
    Base.setindex!(mol::AbstractMolGraph, v::Integer) -> AbstractElement
    Base.setindex!(mol::AbstractMolGraph, e::Edge) -> AbstractElement
    Base.setindex!(mol::AbstractMolGraph, u::Integer, v::Integer) -> AbstractElement

Return properties (vertex or edge attributes).
"""

Base.setindex!(mol::AbstractMolGraph, prop::AbstractElement, v::Integer) = setindex!(mol.vprops, prop, v)
Base.setindex!(mol::AbstractMolGraph, prop::AbstractElement, e::Edge) = setindex!(mol.eprops, prop, e)

# old accessors (deprecated)
props(mol::AbstractMolGraph, v::Integer) = mol[v]
props(mol::AbstractMolGraph, e::Edge) = mol[e]
get_prop(mol::AbstractMolGraph, v::Integer, prop::Symbol) = mol[v][prop]
get_prop(mol::AbstractMolGraph, e::Edge, prop::Symbol) = mol[e][prop]


vpropiter(mol::AbstractMolGraph) = Iterators.map(mol.vprops) do r
    first(r).key, last(r)
end

vpropiter(x::AbstractElement) = Iterators.map(x.vprops) do r
    first(r).key, last(r)
end

epropiter(mol::AbstractMolGraph) = Iterators.map(mol.eprops) do r
    first(r).key, last(r)
end


"""
    SimpleMolGraph{T<:Integer} <: AbstractMolGraph{T}

The base class of molecular graph models based on SimpleGraph
"""
abstract type SimpleMolGraph{T<:Integer} <: AbstractMolGraph{T} end


# SimpleMolGraph interface

Graphs.add_edge!(mol::SimpleMolGraph, u::Integer, v::Integer, prop::AbstractElement
    ) = add_u_edge!(mol, u_edge(mol, u, v), prop)
Graphs.add_edge!(mol::SimpleMolGraph, e::Edge, prop::AbstractElement
    ) = add_edge!(mol, src(e), dst(e), prop)
Graphs.rem_edge!(mol::SimpleMolGraph, u::Integer, v::Integer) = rem_u_edge!(mol, u_edge(mol, u, v))
Graphs.rem_edge!(mol::SimpleMolGraph, e::Edge) = rem_edge!(mol, src(e), dst(e))
Graphs.induced_subgraph(mol::SimpleMolGraph, vlist::AbstractVector{T}
    ) where {T<:Integer} = _induced_subgraph(mol, vlist)
Graphs.induced_subgraph(mol::SimpleMolGraph, elist::AbstractVector{T}
    ) where {T<:Edge} = _induced_subgraph(mol, elist)

function Base.show(io::IO, ::MIME"text/plain", g::SimpleMolGraph)
    print(io, "{$(nv(g)), $(ne(g))} simple molecular graph $(typeof(g))")
end


"""
    u_edge(::Type{T}, src, dst) where T <: Integer -> Edge{T}
    u_edge(g::SimpleGraph{T}, src, dst) where T -> Edge{T}
    u_edge(mol::SimpleMolGraph{T}, src, dst) where T -> Edge{T}

A workaround for UndirectedEdge that are not yet implemented in SimpleGraph
"""
u_edge(::Type{T}, src::T, dst::T) where T<:Integer = src < dst ? Edge{T}(src, dst) : Edge{T}(dst, src)
u_edge(g::SimpleGraph{T}, src::T, dst::T) where T<:Integer = u_edge(T, src, dst)
u_edge(g::SimpleGraph, e::Edge) = u_edge(g, src(e), dst(e))
u_edge(e::Edge{T}) where T<:Integer = u_edge(T, src(e), dst(e))
u_edge(mol::SimpleMolGraph{T}, src::T, dst::T) where T<:Integer = u_edge(T, src, dst)
u_edge(mol::SimpleMolGraph, e::Edge) = u_edge(mol.graph, e)



"""
    edge_neighbors(g::SimpleGraph{T}, u::Integer, v::Integer) where T -> Tuple{Vector{T},Vector{T}}
    edge_neighbors(g::SimpleGraph{T}, e::Edge) where T -> Tuple{Vector{T},Vector{T}}

Return neighbors of the source and destination vertices of the edge, respectively.
"""
edge_neighbors(g::SimpleGraph, u::Integer, v::Integer) = (
    filter(n -> n != v, neighbors(g, u)),
    filter(n -> n != u, neighbors(g, v))
)
edge_neighbors(g::SimpleGraph, e::Edge) = edge_neighbors(g, src(e), dst(e))
edge_neighbors(mol::SimpleMolGraph, u, v) = edge_neighbors(mol.graph, u, v)
edge_neighbors(mol::SimpleMolGraph, e) = edge_neighbors(mol.graph, e)


# if specific index order is required.
# neighbors in SimpleGraph guarantees output index in lexicographic order now,
# but it is not formulated in API
ordered_neighbors = neighbors
ordered_edge_neighbors = edge_neighbors


Base.getindex(mol::SimpleMolGraph{T}, u::T, v::T) where T = mol.eprops[u_edge(T, u, v)]
Base.setindex!(mol::SimpleMolGraph{T}, prop::AbstractElement, u::T, v::T
    ) where T = setindex!(mol.eprops, prop, u_edge(T, u, v))

# old accessors (deprecated)
props(mol::SimpleMolGraph, u::Integer, v::Integer) = mol[u, v]
get_prop(mol::SimpleMolGraph, u::Integer, v::Integer, prop::Symbol) = mol[u, v][prop]
get_prop(mol::SimpleMolGraph, prop::Symbol) = mol[prop]
has_prop(mol::SimpleMolGraph, prop::Symbol) = hasproperty(mol.gprops, prop)
# set_prop! is not available because MolProperty types should be immutable.

# Internally called by descriptor methods (e.g. `is_ring_aromatic(mol)`)
# Do not expose

get_descriptor(mol::SimpleMolGraph, field::Symbol
    ) = getproperty(mol[:descriptors], field)
set_descriptor!(mol::SimpleMolGraph, field::Symbol, value
    ) = setproperty!(mol[:descriptors], field, value)
has_descriptor(mol::SimpleMolGraph, field::Symbol
    ) = hasproperty(mol[:descriptors], field)


"""
    ReactiveMolGraph{T<:Integer,V<:AbstractElement,E<:AbstractElement} <: SimpleMolGraph{T}

The base class of molecule model which have auto-update mechanism of properties.

Typically `ReactiveMolGraph` should have the following properties:
- `graph`: molecular graph topology in `Graphs.SimpleGraph`.
- `vprops`: `Vector` of atom properties (e.g. SDFAtom, SMILESAtom)
- `eprops`: `Vector` of bond properties (e.g. SDFBond, SMILESBond)
- `gprops`: graph-level properties and stored descriptors (e.g. stereocenter)
- `states`: update flags and callback functions for `reactive` property update

"""
abstract type ReactiveMolGraph{T<:Integer,V<:AbstractElement,E<:AbstractElement} <: SimpleMolGraph{T} end

Base.:(==)(g::ReactiveMolGraph, h::ReactiveMolGraph
    ) = g.graph == h.graph && g.vprops == h.vprops && g.eprops == h.eprops && g.gprops == h.gprops

vproptype(::Type{<:ReactiveMolGraph{T,V,E}}) where {T,V,E} = V
eproptype(::Type{<:ReactiveMolGraph{T,V,E}}) where {T,V,E} = E

Base.copy(mol::T) where T <: ReactiveMolGraph = T(
    copy(mol.graph), copy(mol.vprops), copy(mol.eprops), copy(mol.gprops), copy(mol.state))




# Reactions

@kwdef mutable struct ReactionProperty <: AbstractProperty
    # Reaction-level metadata properties (e.g. SDFile option fields)
    metadata::OrderedDict{String,String} = OrderedDict()
    # Parse errors
    logs::Dict{String,String} = Dict()
end


"""
    Reaction{T<:AbstractMolGraph}

Reaction type.
"""
struct Reaction{T<:AbstractMolGraph} <: AbstractReaction
    reactants::Vector{T}
    products::Vector{T}
    rprops::ReactionProperty
end

Reaction{T}() where T = Reaction{T}([], [], ReactionProperty())

Base.eltype(rxn::Reaction) = eltype(rxn.reactants)
Base.eltype(::Type{Reaction{T}}) where T = T