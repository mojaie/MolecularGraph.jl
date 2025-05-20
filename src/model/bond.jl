#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    SDFBond, SMILESBond


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

    function SDFBond(order=1, notation=0, isordered=true)
        new(order, notation, isordered)
    end
end

SDFBond(d::Dict{T,Any}) where T <: Union{AbstractString,Symbol} = SDFBond(
    d[T("order")], d[T("notation")], d[T("isordered")])
SDFBond(arr::Vector) = SDFBond(arr...)

Base.getindex(b::SDFBond, prop::Symbol) = getproperty(b, prop)
Base.:(==)(b1::SDFBond, b2::SDFBond) = all(
    [getfield(b1, f1) == getfield(b2, f2) for (f1, f2) in zip(fieldnames(typeof(b1)), fieldnames(typeof(b2)))])
Base.hash(b::SDFBond, h::UInt
    ) = hash(b.order, hash(b.notation, hash(b.isordered, h)))

to_dict(::Val{:standard}, b::SDFBond) = Any[b.order, b.notation, b.isordered]


"""
    SMILESBond

SMILES bond property type.
"""
struct SMILESBond
    order::Int
    isaromatic::Bool
    direction::Symbol  # :up, :down or :unspecified

    function SMILESBond(order=1, isaromatic=false, direction=:unspecified)
        new(order, isaromatic, direction)
    end
end

SMILESBond(d::Dict{T,Any}) where T <: Union{AbstractString,Symbol} = SMILESBond(
    d[T("order")], d[T("isaromatic")], Symbol(d[T("direction")]))
SMILESBond(arr::Vector) = SMILESBond(arr[1], arr[2], Symbol(arr[3]))

Base.getindex(b::SMILESBond, prop::Symbol) = getproperty(b, prop)
Base.:(==)(b1::SMILESBond, b2::SMILESBond) = all(
    [getfield(b1, f1) == getfield(b2, f2) for (f1, f2) in zip(fieldnames(typeof(b1)), fieldnames(typeof(b2)))])
Base.hash(b::SMILESBond, h::UInt
    ) = hash(b.order, hash(b.isaromatic, hash(b.direction, h)))

to_dict(::Val{:standard}, b::SMILESBond) = Any[b.order, b.isaromatic, b.direction]
