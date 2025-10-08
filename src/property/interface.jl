#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#


# Common interfaces of AbstractProperty


# Metadata

reconstruct(::Val{:metadata}, ::Type{T}, data::JSON.Object{String,Any}
    ) where T <: AbstractProperty = OrderedDict(d[1] => d[2] for d in data)

to_dict(::Val{:metadata}, ::Val{:default}, gprop::AbstractProperty
    ) = Any[collect(d) for d in gprop.metadata]


# Metadata specific shorthands (e.g. mol["compound_id"] = "CP000001")

Base.getindex(mol::SimpleMolGraph, key::String) = mol.gprops.metadata[key]

function Base.setindex!(mol::SimpleMolGraph, value::String, key::String)
    # skip dispatch (Metadata update would not affect graph state)
    mol.gprops.metadata[key] = value
end

# old accessors (deprecated)
get_prop(mol::SimpleMolGraph, key::String) = mol[key]
has_prop(mol::SimpleMolGraph, key::String) = haskey(mol.gprops.metadata, key)
set_prop!(mol::SimpleMolGraph, key::String, value::String) = setindex!(mol, value, key)


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

Base.copy(desc::T) where T <: MolDescriptor = T(
    copy_vec_of_vec(desc.sssr), copy(desc.lone_pair), copy(desc.apparent_valence),
    copy(desc.valence), copy(desc.is_ring_aromatic), copy(desc.atom_charge),
    copy(desc.bond_order), copy_vec_of_vec(desc.coords2d), copy_vec_of_vec(desc.coords3d),
    copy_vec_of_vec(desc.draw2d_bond_style)
)

# `reconstruct` should be interfaced indivisually

to_dict(::Val{:descriptors}, ::Val{:default}, gprop::AbstractProperty
    ) = to_dict(Val(:default), gprop.descriptors)


# Properties

"""
    MolProperty{T} <: SimpleMolProperty{T}

Container of graph-level molecule properties compatible with `ReactiveMolGraph`.
"""
@kwdef struct MolProperty{T} <: SimpleMolProperty{T}
    stereocenter::StereocenterMap{T} = StereocenterMap{T}()
    stereobond::StereobondMap{T} = StereobondMap{T}()
    pyrrole_like::Vector{T} = T[]  # pyrrole H position for SMILES kekulization
    smarts_input::String = ""
    smarts_lexical_succ::Vector{Vector{T}} = Vector{T}[]  # lexical index used for stereochem
    descriptors::MolDescriptor{T} = MolDescriptor{T}()
    # Graph-level metadata properties (e.g. SDFile option fields)
    metadata::OrderedDict{String,String} = OrderedDict{String,String}()
    # Parse errors
    logs::Dict{String,String} = Dict{String,String}()
end


Base.copy(prop::T) where T <: MolProperty = T(
    copy(prop.stereocenter), copy(prop.stereobond), copy(prop.pyrrole_like),
    prop.smarts_input, copy_vec_of_vec(prop.smarts_lexical_succ), copy(prop.descriptors),
    copy(prop.metadata), copy(prop.logs)
)

reconstruct(::Val{:pyrrole_like}, ::Type{MolProperty{T}}, data::Vector{Any}) where T = T.(data)
reconstruct(::Val{:smarts_input}, ::Type{MolProperty{T}}, data::String) where T = data
reconstruct(::Val{:smarts_lexical_succ}, ::Type{MolProperty{T}}, data::Vector{Any}
    ) where T = [T.(d) for d in data]
reconstruct(::Val{:descriptors}, ::Type{MolProperty{T}}, data::JSON.Object{String,Any}
    ) where T = reconstruct(MolDescriptor{T}, data)



