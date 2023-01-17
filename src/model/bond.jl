#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    SDFBond, SMILESBond,
    setorder, setnotation, setstereo


struct SDFBond
    """Bond
    * Notation
        * Single bond
            * 0: u - v
            * 1: u ◀ v (Up-arrow)
            * 2: v ◀ u -> for SMILES coordgen compatibility # deprecated
            * 4: u ~ v (Up or down)
            * 6: u ◁ v (Down-arrow)
            * 7: v ◁ u -> for SMILES coordgen compatibility # deprecated
        * Double bond
            * 0: v ニ u (clockwise, default)
            * 1: u ニ v (counter-clockwise) # deprecated
            * 2: u ＝ v (equal length, for terminal bond by default) # deprecated
            * 3: u × v (Cis-Trans Unknown)
    """
    order::Int
    notation::Int
    is_ordered::Bool
    stereo::Symbol  # deprecated
end

SDFBond(order, notation, is_ordered) = SDFBond(order, notation, is_ordered, :unspecified)
SDFBond() = SDFBond(1, 0, true, :unspecified)
SDFBond(d::Dict{T,Any}) where T <: Union{AbstractString,Symbol} = SDFBond(
    d[T("order")], d[T("notation")], d[T("is_ordered")])

Base.getindex(b::SDFBond, prop::Symbol) = getproperty(b, prop)

function to_dict(b::SDFBond)
    data = Dict{String,Any}()
    for field in fieldnames(SDFBond)
        data[string(field)] = getfield(b, field)
    end
    return data
end

setorder(b::SDFBond, order
    ) = SDFBond(order, b.notation, b.is_ordered, b.stereo)  # deprecated
setnotation(b::SDFBond, notation
    ) = SDFBond(b.order, notation, b.is_ordered, b.stereo)  # deprecated
setstereo(edge::SDFBond, cistrans
    ) = SDFBond(b.order, b.notation, b.is_ordered, cistrans)  # deprecated


struct SMILESBond
    order::Int
    is_aromatic::Bool
    direction::Symbol  # :up or :down
    stereo::Symbol  # deprecated
end

SMILESBond(order, is_aromatic, direction) = SMILESBond(order, is_aromatic, direction, :unspecified)
SMILESBond() = SMILESBond(1, false, :unspecified)
SMILESBond(d::Dict{T,Any}) where T <: Union{AbstractString,Symbol} = SDFBond(
    d[T("order")], d[T("is_aromatic")], d[T("direction")])

Base.getindex(b::SMILESBond, prop::Symbol) = getproperty(b, prop)

function todict(b::SMILESBond)
    data = Dict{String,Any}()
    for field in fieldnames(SMILESBond)
        data[string(field)] = getfield(b, field)
    end
    return data
end

setorder(b::SMILESBond, order
    ) = SMILESBond(order, b.is_aromatic, b.direction, b.stereo)  # deprecated
setstereo(b::SMILESBond, cistrans
    ) = SMILESBond(b.order, b.is_aromatic, b.direction, cistrans)  # deprecated
