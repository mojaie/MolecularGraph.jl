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
mutable struct MolDescriptor{T} <: SimpleMolProperty{T}
    # standardized atom charges and bond orders
    atom_charge::Vector{Int}
    bond_order::Vector{Int}
    # cached descriptors that are expensive or frequently used
    sssr::Vector{Vector{T}}
    lone_pair::Vector{Int}
    apparent_valence::Vector{Int}
    valence::Vector{Int}
    is_ring_aromatic::Vector{Bool}
    # coordinates
    coords2d::Vector{Coords2d}
    coords3d::Vector{Coords3d}
    draw2d_bond_style::Vector{Vector{Symbol}}  # wedge notation in drawing
end

function MolDescriptor{T}(;
        atom_charge::Vector{Int}=Int[],
        bond_order::Vector{Int}=Int[],
        sssr::Vector{Vector{T}}=Vector{T}[],
        lone_pair::Vector{Int}=Int[],
        apparent_valence::Vector{Int}=Int[],
        valence::Vector{Int}=Int[],
        is_ring_aromatic::Vector{Bool}=Bool[],
        coords2d::Vector{Coords2d}=Coords2d[],
        coords3d::Vector{Coords3d}=Coords3d[],
        draw2d_bond_style::Vector{Vector{Symbol}}=Vector{Symbol}[]) where T <: Integer
    return MolDescriptor{T}(
        atom_charge, bond_order, sssr, lone_pair, apparent_valence,
        valence, is_ring_aromatic, coords2d, coords3d, draw2d_bond_style
    )
end

Base.copy(desc::T) where T <: MolDescriptor = T(
    copy(desc.atom_charge), copy(desc.bond_order), copy_vec_of_vec(desc.sssr),
    copy(desc.lone_pair), copy(desc.apparent_valence),
    copy(desc.valence), copy(desc.is_ring_aromatic),
    [copy(c) for c in desc.coords2d], [copy(c) for c in desc.coords3d],
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
struct MolProperty{T} <: SimpleMolProperty{T}
    stereocenter::StereocenterMap{T}
    stereobond::StereobondMap{T}
    pyrrole_like::Vector{T}  # pyrrole H position for SMILES kekulization
    smarts_input::String  # TODO: should be metadata
    smarts_lexical_succ::Vector{Vector{T}}  # lexical index used for stereochem
    descriptors::MolDescriptor{T}
    metadata::OrderedDict{String,String}  # e.g. SDFile option fields
    logs::Dict{String,String}  # Parse errors
end

function MolProperty{T}(;
        stereocenter::StereocenterMap{T}=StereocenterMap{T}(),
        stereobond::StereobondMap{T}=StereobondMap{T}(),
        pyrrole_like::Vector{T}=T[],
        smarts_input::String="",
        smarts_lexical_succ::Vector{Vector{T}}=Vector{T}[],
        descriptors::MolDescriptor{T}=MolDescriptor{T}(),
        metadata::OrderedDict{String,String}=OrderedDict{String,String}(),
        logs::Dict{String,String}=Dict{String,String}()) where T <: Integer
    return MolProperty{T}(
        stereocenter, stereobond, pyrrole_like, smarts_input, smarts_lexical_succ,
        descriptors, metadata, logs
    )
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



