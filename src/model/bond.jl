#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    SDFBond, SMILESBond


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
    isordered::Bool
    stereo::Symbol  # deprecated

    function SDFBond(order=1, notation=0, isordered=true)
        new(order, notation, isordered, :unspecified)
    end
end

SDFBond(d::Dict{T,Any}) where T <: Union{AbstractString,Symbol} = SDFBond(
    d[T("order")], d[T("notation")], d[T("isordered")])

Base.getindex(b::SDFBond, prop::Symbol) = getproperty(b, prop)

function to_dict(b::SDFBond)
    data = Dict{String,Any}()
    for field in fieldnames(SDFBond)
        data[string(field)] = getfield(b, field)
    end
    return data
end



struct SMILESBond
    order::Int
    isaromatic::Bool
    direction::Symbol  # :up or :down
    stereo::Symbol  # deprecated

    function SMILESBond(order=1, isaromatic=false, direction=:unspecified)
        new(order, isaromatic, direction, :unspecified)
    end
end

SMILESBond(d::Dict{T,Any}) where T <: Union{AbstractString,Symbol} = SMILESBond(
    get(d, T("order"), 1),
    get(d, T("isaromatic"), false),
    get(d, T("direction"), :unspecified)
)

Base.getindex(b::SMILESBond, prop::Symbol) = getproperty(b, prop)

function todict(b::SMILESBond)
    data = Dict{String,Any}()
    for field in fieldnames(SMILESBond)
        data[string(field)] = getfield(b, field)
    end
    return data
end
