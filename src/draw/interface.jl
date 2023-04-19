#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

using Printf

export
    Canvas,
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
    atom_annotation!,

    Color,
    DEFAULT_ATOM_COLOR,
    RASMOL_ATOM_COLOR,

    SvgCanvas,
    tosvg,
    drawsvg,
    initcanvas!,

    html_fixed_size,
    html_grid,

    atom_color, is_atom_visible,
    single_bond_style, double_bond_style, bond_style,
    chargesign, atommarkup, atomhtml,
    draw2d!, drawatomindex!, sethighlight!

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
        color::Color, implicithconnected::Int, charge::Int
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
