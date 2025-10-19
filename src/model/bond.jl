#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

# Bond interfaces


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
struct SDFBond <: StandardBond
    order::Int
    notation::Int
    isordered::Bool

    function SDFBond(
            order::Int,
            notation::Int=0,
            isordered::Bool=true)
        order > 3 && error("sdfile parse error - unsupported bond order $(order)")
        new(order, notation, isordered)
    end
end

function SDFBond(;
        order::Int=1,
        notation::Int=0,
        isordered::Bool=true)
    return SDFBond(order, notation, isordered)
end

SDFBond(d::Dict{Symbol,Any}) = SDFBond(; NamedTuple((k, v) for (k, v) in d)...)

bond_order(b::SDFBond) = b.order

StructUtils.structlike(::StructUtils.StructStyle, ::Type{SDFBond}) = false

function JSON.lower(x::SDFBond)
    rcd = Dict{String,Any}()
    x.order == 1 || setindex!(rcd, x.order, "order")
    x.notation == 0 || setindex!(rcd, x.notation, "notation")
    x.isordered === true || setindex!(rcd, x.isordered, "isordered")
    return rcd
end

JSON.lift(::Type{SDFBond}, x) = SDFBond(; NamedTuple((Symbol(k), v) for (k, v) in x)...)



"""
    SMILESBond

SMILES bond property type.
"""
struct SMILESBond <: StandardBond
    order::Int
    isaromatic::Bool
    direction::Symbol  # :up, :down or :unspecified

    function SMILESBond(
            order::Int,
            isaromatic::Bool=false,
            direction::Union{AbstractString,Symbol}=:unspecified)
        new(order, isaromatic, Symbol(direction))
    end
end

function SMILESBond(;
        order::Int=1,
        isaromatic::Bool=false,
        direction::Union{AbstractString,Symbol}=:unspecified)
    return SMILESBond(order, isaromatic, Symbol(direction))
end

SMILESBond(d::Dict{Symbol,Any}) = SMILESBond(; NamedTuple((k, v) for (k, v) in d)...)

bond_order(b::SMILESBond) = b.order

StructUtils.structlike(::StructUtils.StructStyle, ::Type{SMILESBond}) = false

function JSON.lower(x::SMILESBond)
    rcd = Dict{String,Any}()
    x.order == 1 || setindex!(rcd, x.order, "order")
    x.isaromatic === false || setindex!(rcd, x.isaromatic, "isaromatic")
    x.direction === :unspecified || setindex!(rcd, x.direction, "direction")
    return rcd
end

JSON.lift(::Type{SMILESBond}, x) = SMILESBond(; NamedTuple((Symbol(k), v) for (k, v) in x)...)



"""
    CommonChemBond

CommonChem bond property type.
"""
struct CommonChemBond <: StandardBond
    type::Int  # bond order

    function CommonChemBond(
            type::Union{Int,Nothing},
            bo::Union{Int,Nothing}=nothing)
        # RDkit notation is 'bo', not 'type'
        new(something(bo, type))
    end
end

function CommonChemBond(;
        type::Union{Int,Nothing}=1,
        bo::Union{Int,Nothing}=nothing)
    # RDkit notation is 'bo', not 'type'
    CommonChemBond(something(bo, type))
end

CommonChemBond(d::Dict{Symbol,Any}) = CommonChemBond(; NamedTuple((k, v) for (k, v) in d)...)
CommonChemBond(x::SDFBond) = CommonChemBond(x.order)
CommonChemBond(x::SMILESBond) = CommonChemBond(x.order)

bond_order(b::CommonChemBond) = b.type

StructUtils.structlike(::StructUtils.StructStyle, ::Type{CommonChemBond}) = false

function JSON.lower(x::CommonChemBond)
    rcd = Dict{String,Any}()
    x.type == 1 || setindex!(rcd, x.type, "bo")
    return rcd
end

JSON.lift(::Type{CommonChemBond}, x) = CommonChemBond(; NamedTuple((Symbol(k), v) for (k, v) in x)...)
