#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@kwdef mutable struct Descriptors{T}
    # cached relatively expensive descriptors
    sssr::Vector{Vector{T}} = Vector{T}[]
    lone_pair::Vector{Int} = Int[]
    apparent_valence::Vector{Int} = Int[]
    valence::Vector{Int} = Int[]
    is_ring_aromatic::Vector{Bool} = Bool[]
    # standardized atom charges and bond orders
    atom_charge::Vector{Int} = Int[]
    bond_order::Vector{Int} = Int[]
end

function Descriptors{T}(data::Dict{String,Any}) where T
    desc = Descriptors{T}()
    for (k, v) in data
        setproperty!(desc, k, v)
    end
    return desc
end

function Base.:(==)(g::Descriptors, h::Descriptors)
    for sym in fieldnames(typeof(g))
        getproperty(g, sym) == getproperty(h, sym) || return false
    end
    return true
end



@kwdef mutable struct MolGraphProperty{T}
    stereocenter::Dict{T,Tuple{T,T,T,Bool}} = Dict{T,Tuple{T,T,T,Bool}}()
    stereobond::Dict{Edge{T},Tuple{T,T,Bool}} = Dict{Edge{T},Tuple{T,T,Bool}}()
    pyrrole_like::Vector{T} = T[]  # to keep pyrrole H position
    smarts_input::String = ""
    smarts_lexical_succ::Vector{Vector{T}} = Vector{T}[]  # lexical index used for stereochem
    smarts_connectivity::Vector{Vector{T}} = Vector{T}[]  # SMARTS connectivity query
    descriptors::Descriptors{T} = Descriptors{T}()
    coords2d::Vector{Vector{Point2d}} = Vector{Point2d}[]
    draw2d_bond_style::Vector{Vector{Symbol}} = Vector{Symbol}[]  # wedge notation
    coords3d::Vector{Vector{Point3d}} = Vector{Point3d}[]
    # Graph-level metadata properties (e.g. SDFile option fields)
    metadata::OrderedDict{String,String} = OrderedDict{String,String}()
    # Parse errors
    logs::Dict{String,String} = Dict{String,String}()
end

function MolGraphProperty{T}(data::Dict{String,Any}) where T
    gprop = MolGraphProperty{T}()
    for sym in fieldnames(typeof(gprop))
        reconstruct!(Val(sym), gprop, data[string(sym)])
    end
    return gprop
end

function Base.:(==)(g::MolGraphProperty, h::MolGraphProperty)
    for sym in fieldnames(typeof(g))
        getproperty(g, sym) == getproperty(h, sym) || return false
    end
    return true
end

function to_dict(::Val{T}, gprop::MolGraphProperty) where T
    data = Dict{String,Any}()
    for sym in fieldnames(typeof(gprop))
        data[string(sym)] = to_dict(Val(T), Val(sym), gprop)
    end
    return data
end


function reconstruct!(::Val{T}, gprop::MolGraphProperty, data) where T
    setproperty!(gprop, T, reconstruct(Val(T), gprop, data))
end


function remap!(::Val{T}, gprop::MolGraphProperty, vmap::Dict) where T
    # vmap[old] -> new
    setproperty!(gprop, T, remap(Val(T), gprop, vmap))
end


function remap_gprops(mol::ReactiveMolGraph{T,V,E}, vmap::Dict{T,T}) where {T,V,E}
    gprop = MolGraphProperty{T}()
    for k in fieldnames(typeof(mol.gprops))
        setproperty!(gprop, k, remap(Val(k), mol.gprops, vmap))
    end
    return gprop
end

function remap_gprops!(mol::ReactiveMolGraph, vmap::Dict)
    for k in fieldnames(typeof(mol.gprops))
        remap!(Val(k), mol.gprops, vmap)
    end
end


to_dict(
    ::Val{:default}, ::Val{T}, gprop::MolGraphProperty) where T = getproperty(gprop, T)
reconstruct(::Val{T}, gprop::MolGraphProperty, data) where T = data
remap(::Val{T}, gprop::MolGraphProperty, vmap::Dict) where T = getproperty(gprop, T)


function to_dict(::Val{:default}, ::Val{:descriptors}, gprop::MolGraphProperty)
    data = Dict{String,Any}()
    for sym in fieldnames(typeof(gprop.descriptors))
        data[string(sym)] = getproperty(gprop.descriptors, sym)
    end
    return data
end

function reconstruct(::Val{:descriptors}, gprop::MolGraphProperty{T}, data) where T
    desc = Descriptors{T}()
    for sym in fieldnames(typeof(gprop.descriptors))
        setproperty!(desc, sym, data[string(sym)])
    end
    return desc
end

# Descriptors would be recalculated by the update callback. Do nothing.
remap(::Val{:descriptors}, gprop::MolGraphProperty, vmap::Dict) = gprop.descriptors


to_dict(
    ::Val{:default}, ::Val{:metadata}, gprop::MolGraphProperty) = [collect(d) for d in gprop.metadata]
reconstruct(::Val{:metadata}, gprop::MolGraphProperty, data) = OrderedDict(d[1] => d[2] for d in data)
remap(::Val{:metadata}, gprop::MolGraphProperty, vmap::Dict) = gprop.metadata


# Metadata shortcuts

function set_prop!(mol::ReactiveMolGraph, prop::String, value::String)
    # Metadata update would not affect graph state
    mol.gprops.metadata[prop] = value
end

get_prop(mol::ReactiveMolGraph, prop::String) = mol.gprops.metadata[prop]
has_prop(mol::ReactiveMolGraph, prop::String) = haskey(mol.gprops.metadata, prop)

Base.getindex(mol::ReactiveMolGraph, key::String) = get_prop(mol, key)
Base.setindex!(mol::ReactiveMolGraph, value::String, key::String) = set_prop!(mol, key, value)