#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#


# Common interfaces of AbstractProperty

"""
    remap!(container::SimpleMolProperty{T}, vmap::Dict{T,T}) where T <: Integer

Remap vertices according to the given vmap (new_v = old_v -> vmap[old_v]).
"""
function remap!(container::SimpleMolProperty{T}, vmap::Dict{T,T}) where T <: Integer
    # vmap[old] -> new
    for sym in fieldnames(typeof(container))
        remap!(Val(sym), container, vmap)
    end
    return
end

function remap!(::Val, container::SimpleMolProperty{T}, vmap::Dict{T,T}) where T <: Integer
    return
end


"""
    reconstruct(::Type{T}, data::Dict{String,Any}) where T <: AbstractProperty

Reconstruct properties from a dict, typically form a serialized form of the data (e.g. JSON)
"""
function reconstruct(::Type{T}, data::Dict{String,Any}) where T <: AbstractProperty
    container = T()
    for sym in fieldnames(typeof(container))
        setproperty!(container, sym, reconstruct(Val(sym), T, data[string(sym)]))
    end
    return container
end

reconstruct(::Val{V}, ::Type{T}, @nospecialize(data)) where {V,T<:AbstractProperty} = data


"""
    reconstruct(::Type{T}, data::Dict{String,Any}) where T <: AbstractProperty

Dump properties to a dict for serialization.
"""
function to_dict(::Val{F}, container::AbstractProperty) where F
    data = Dict{String,Any}()
    for sym in fieldnames(typeof(container))
        data[string(sym)] = to_dict(Val(sym), Val(F), container)
    end
    return data
end

to_dict(::Val{V}, ::Val{F}, container::AbstractProperty) where {V,F} = getproperty(container, V)


function Base.:(==)(g::T, h::T) where T <: AbstractProperty
    for sym in fieldnames(typeof(g))
        getproperty(g, sym) == getproperty(h, sym) || return false
    end
    return true
end


# Metadata

reconstruct(::Val{:metadata}, ::Type{T}, @nospecialize(data)
    ) where T <: AbstractProperty = OrderedDict(d[1] => d[2] for d in data)

to_dict(::Val{:metadata}, ::Val{:default}, gprop::AbstractProperty
    ) = Any[collect(d) for d in gprop.metadata]


# Metadata specific shorthands (e.g. mol["compound_id"] = "CP000001")

function set_prop!(mol::ReactiveMolGraph, prop::String, value::String)
    # Metadata update would not affect graph state
    mol.gprops.metadata[prop] = value
end

get_prop(mol::ReactiveMolGraph, prop::String) = mol.gprops.metadata[prop]
has_prop(mol::ReactiveMolGraph, prop::String) = haskey(mol.gprops.metadata, prop)

Base.getindex(mol::ReactiveMolGraph, key::String) = get_prop(mol, key)
Base.setindex!(mol::ReactiveMolGraph, value::String, key::String) = set_prop!(mol, key, value)


# Descriptors

"""
    MolDescriptor{T} <: SimpleMolProperty{T}

Container of calculated secondary properties (descriptor) compatible with `ReactiveMolGraph`.
"""
@kwdef mutable struct MolDescriptor{T} <: SimpleMolProperty{T}
    # cached relatively expensive descriptors
    sssr::Vector{Vector{T}} = Vector{T}[]
    lone_pair::Vector{Int} = Int[]
    apparent_valence::Vector{Int} = Int[]
    valence::Vector{Int} = Int[]
    is_ring_aromatic::Vector{Bool} = Bool[]
    # standardized atom charges and bond orders
    atom_charge::Vector{Int} = Int[]
    bond_order::Vector{Int} = Int[]
    # coordinates
    coords2d::Vector{Vector{Point2d}} = Vector{Point2d}[]
    coords3d::Vector{Vector{Point3d}} = Vector{Point3d}[]
    # wedge notation in drawing
    draw2d_bond_style::Vector{Vector{Symbol}} = Vector{Symbol}[]
end

# `reconstruct` should be interfaced indivisually

to_dict(::Val{:descriptors}, ::Val{:default}, gprop::AbstractProperty
    ) = to_dict(Val(:default), gprop.descriptors)

function remap!(::Val{:descriptors}, gprop::SimpleMolProperty{T},
        vmap::Dict{T,T}) where T <: Integer
    for sym in fieldnames(typeof(gprop.descriptors))
        remap!(Val(sym), gprop.descriptors, vmap)
    end
    return
end


# Properties

"""
    MolProperty{T} <: SimpleMolProperty{T}

Container of graph-level molecule properties compatible with `ReactiveMolGraph`.
"""
@kwdef mutable struct MolProperty{T} <: SimpleMolProperty{T}
    stereocenter::Dict{T,Tuple{T,T,T,Bool}} = Dict{T,Tuple{T,T,T,Bool}}()
    stereobond::Dict{Edge{T},Tuple{T,T,Bool}} = Dict{Edge{T},Tuple{T,T,Bool}}()
    pyrrole_like::Vector{T} = T[]  # pyrrole H position for SMILES kekulization
    smarts_input::String = ""
    smarts_lexical_succ::Vector{Vector{T}} = Vector{T}[]  # lexical index used for stereochem
    descriptors::MolDescriptor{T} = MolDescriptor{T}()
    # Graph-level metadata properties (e.g. SDFile option fields)
    metadata::OrderedDict{String,String} = OrderedDict{String,String}()
    # Parse errors
    logs::Dict{String,String} = Dict{String,String}()
end


reconstruct(::Val{:descriptors}, ::Type{MolProperty{T}}, @nospecialize(data)
    ) where T = reconstruct(MolDescriptor{T}, data)



# Interfaces for ReactiveMolGraph

"""
    remap_gprops(mol::ReactiveMolGraph, vmap::Dict{T,T}) -> MolProperty

Return updated `MolProperty`.

Only `Graph.rem_vertex!` variants should call this function to remap vertices
after removal.
"""
function remap_gprops(mol::ReactiveMolGraph, vmap::Dict{T,T}) where T <: Integer
    gprop = deepcopy(mol.gprops)
    remap!(gprop, vmap)
    return gprop
end

"""
    remap_gprops!(mol::ReactiveMolGraph, vmap::Dict{T,T}) -> Nothing

Update `MolProperty` in place.

Only `Graph.rem_vertex!` variants should call this function to remap vertices
after removal.
"""
function remap_gprops!(mol::ReactiveMolGraph, vmap::Dict{T,T}) where T <: Integer
    remap!(mol.gprops, vmap)
    return
end


get_prop(mol::ReactiveMolGraph, prop::Symbol) = getproperty(mol.gprops, prop)
has_prop(mol::ReactiveMolGraph, prop::Symbol) = hasproperty(mol.gprops, prop)
# Note: editing graph-level properties may break consistency.
set_prop!(mol::ReactiveMolGraph, prop::Symbol, value
    ) = setproperty!(mol.gprops, prop, value)

get_descriptor(mol::ReactiveMolGraph, field::Symbol
    ) = getproperty(mol.gprops.descriptors, field)
has_descriptor(mol::ReactiveMolGraph, field::Symbol
    ) = hasproperty(mol.gprops.descriptors, field)
set_descriptor!(mol::ReactiveMolGraph, field::Symbol, value
    ) = setproperty!(mol.gprops.descriptors, field, value)

