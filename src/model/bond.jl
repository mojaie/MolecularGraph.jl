#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#


bond_order(b::AbstractDict) = b[:order]


"""
    SDFBond

SDFile (CTAB) bond property type.

* SDFile bond notation
    * Single bond
        * 0: u - v
        * 1: u ◀ v (Up-arrow)
        * 4: u ~ v (Up or down)
        * 6: u ◁ v (Down-arrow)
    * Double bond
        * 0: v = u
        * 3: u x v (Cis-Trans Unknown)
"""
struct SDFBond
    order::Int
    notation::Int
    isordered::Bool
end

function SDFBond(;
        order::Int=1,
        notation::Int=0,
        isordered::Bool=true)
    return SDFBond(order, notation, isordered)
end

SDFBond(d::Dict{String,Any}
    ) = SDFBond(; NamedTuple((Symbol(k), v) for (k, v) in d)...)
SDFBond(d::Dict{Symbol,Any}
    ) = SDFBond(; NamedTuple((k, v) for (k, v) in d)...)

Base.getindex(b::SDFBond, prop::Symbol) = getproperty(b, prop)
Base.:(==)(b1::SDFBond, b2::SDFBond) = all(
    [getfield(b1, f1) == getfield(b2, f2) for (f1, f2) in zip(fieldnames(typeof(b1)), fieldnames(typeof(b2)))])
Base.hash(b::SDFBond, h::UInt
    ) = hash(b.order, hash(b.notation, hash(b.isordered, h)))

bond_order(b::SDFBond) = b.order

ELEMENT_TYPE_REGISTRY["SDFBond"] = SDFBond

function to_dict(::Val{:default}, b::SDFBond)
    rcd = Dict{String,Any}()
    b.order == 1 || setindex!(rcd, b.order, "order")
    b.notation == 0 || setindex!(rcd, b.notation, "notation")
    b.isordered === true || setindex!(rcd, b.isordered, "isordered")
    return rcd
end

function to_dict(::Val{:rdkit}, b::SDFBond)
    rcd = Dict{String,Any}()
    b.order == 1 || setindex!(rcd, b.order, "bo")
    return rcd
end


"""
    SMILESBond

SMILES bond property type.
"""
struct SMILESBond
    order::Int
    isaromatic::Bool
    direction::Symbol  # :up, :down or :unspecified
end

function SMILESBond(;
        order::Int=1,
        isaromatic::Bool=false,
        direction::Union{AbstractString,Symbol}=:unspecified)
    return SMILESBond(order, isaromatic, Symbol(direction))
end

SMILESBond(d::Dict{String,Any}
    ) = SMILESBond(; NamedTuple((Symbol(k), v) for (k, v) in d)...)
SMILESBond(d::Dict{Symbol,Any}
    ) = SMILESBond(; NamedTuple((k, v) for (k, v) in d)...)

Base.getindex(b::SMILESBond, prop::Symbol) = getproperty(b, prop)
Base.:(==)(b1::SMILESBond, b2::SMILESBond) = all(
    [getfield(b1, f1) == getfield(b2, f2) for (f1, f2) in zip(fieldnames(typeof(b1)), fieldnames(typeof(b2)))])
Base.hash(b::SMILESBond, h::UInt
    ) = hash(b.order, hash(b.isaromatic, hash(b.direction, h)))

bond_order(b::SMILESBond) = b.order

ELEMENT_TYPE_REGISTRY["SMILESBond"] = SMILESBond

function to_dict(::Val{:default}, b::SMILESBond)
    rcd = Dict{String,Any}()
    b.order == 1 || setindex!(rcd, b.order, "order")
    b.isaromatic === false || setindex!(rcd, b.isaromatic, "isaromatic")
    b.direction === :unspecified || setindex!(rcd, b.direction, "direction")
    return rcd
end

function to_dict(::Val{:rdkit}, b::SMILESBond)
    rcd = Dict{String,Any}()
    b.order == 1 || setindex!(rcd, b.order, "bo")
    return rcd
end


"""
    CommonChemBond

CommonChem bond property type.
"""
struct CommonChemBond
    type::Int  # bond order
end

function CommonChemBond(;
        type::Union{Int,Nothing}=1,
        bo::Union{Int,Nothing}=nothing)
    # RDkit notation is 'bo', not 'type'
    CommonChemBond(something(bo, type))
end

CommonChemBond(d::Dict{String,Any}
    ) = CommonChemBond(; NamedTuple((Symbol(k), v) for (k, v) in d)...)
CommonChemBond(d::Dict{Symbol,Any}
    ) = CommonChemBond(; NamedTuple((k, v) for (k, v) in d)...)

CommonChemBond(arr::Vector) = CommonChemBond(arr[1])

Base.getindex(b::CommonChemBond, prop::Symbol) = getproperty(b, prop)
Base.:(==)(b1::CommonChemBond, b2::CommonChemBond) = all(
    [getfield(b1, f1) == getfield(b2, f2) for (f1, f2) in zip(fieldnames(typeof(b1)), fieldnames(typeof(b2)))])
Base.hash(b::CommonChemBond, h::UInt) = hash(b.type, h)

bond_order(b::CommonChemBond) = b.type

ELEMENT_TYPE_REGISTRY["CommonChemBond"] = CommonChemBond

function to_dict(::Val, b::CommonChemBond)
    rcd = Dict{String,Any}()
    b.type == 1 || setindex!(rcd, b.type, "bo")
    return rcd
end