#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    Canvas, Color,
    singlebond!,
    wedged!,
    dashedwedged!,
    wavesingle!,
    doublebond!,
    clockwisedouble!,
    counterdouble!,
    crossdouble!,
    triplebond!,
    atomsymbolright!,
    atomsymbolcenter!,
    atomsymbolleft!,
    atom_annotation!


using MolecularGraph.MolecularGraphGeometry: Segment2D, _point, _vector, _u, _v


struct Color
    r::Int
    g::Int
    b::Int
end

Formatting.format(expr::String, c::Color) = format(expr, c.r, c.g, c.b)


abstract type Canvas end


"""
    setbond!(
        canvas::Canvas, bondorder::Int, bondnotation::Int, coords::Segment2D,
        ucolor::Color, vcolor::Color, u_visible::Bool, v_visible::Bool
    )

Interface for bond drawing.
"""
function setbond! end


"""
    setatomright!(
        canvas::Canvas, coords::Array{Float64,2}, atomsymbol::Symbol,
        color::Color, implicithcount::Int, charge::Int
    )

Interface for atom drawing.
"""
function setatomright! end
function setatomcenter! end
function setatomleft! end


"""
    setatomnote!(
        canvas::Canvas, coords::Array{Float64,2}, text::String, color::Color,
        bgcolor::Color
    )

Interface for atom note drawing.
"""
function setatomnote! end
