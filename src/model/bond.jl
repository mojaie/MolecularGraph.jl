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

Base.getindex(b::SDFBond, prop::Symbol) = getproperty(b, prop)
to_dict(b::SDFBond) = Dict{String,Any}(
    string(field) => getfield(b, field) for field in fieldnames(SDFBond))


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
    get(d, T("order"), 1),
    get(d, T("isaromatic"), false),
    get(d, T("direction"), :unspecified)
)

Base.getindex(b::SMILESBond, prop::Symbol) = getproperty(b, prop)
to_dict(b::SMILESBond) = Dict{String,Any}(
    string(field) => getfield(b, field) for field in fieldnames(SMILESBond))

