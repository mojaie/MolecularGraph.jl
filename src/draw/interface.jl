#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    drawsvg,
    drawpng,
    html_fixed_size,
    html_grid


abstract type Canvas end


# Custom 2D drawing systems using `draw2d!` should implement following methods (see ./svg.jl and ./cairo.jl)

"""
    initcanvas!(canvas::Canvas, coords::AbstractArray{Float64}, boundary::Tuple)

Convert the coordinate system of SDFile to that of the graphics library and set the drawing area.
"""
function initcanvas! end


"""
    atommarkupleft(
        canvas::Canvas, atomsymbol::Symbol, charge::Int, implicith::Int) -> String

Return a HTML or SVG text of the atom and implicit hydrogens in a right-to-left direction.

The number of hydrogens is expressed as a subscript and the number of charges as a superscript.
"""
function atommarkupleft end


"""
    atommarkupright(
        canvas::Canvas, atomsymbol::Symbol, charge::Int, implicith::Int) -> String

Return a HTML or SVG text of the atom and implicit hydrogens in a left-to-right direction.

The number of hydrogens is expressed as a subscript and the number of charges as a superscript.
"""
function atommarkupright end


"""
    drawtextleft!(
        canvas::Canvas, coords::Array{Float64,2}, text::String, color::Color)

Draw characters of an atom and implicit hydrogens in a right-to-left direction.
"""
function drawtextleft! end


"""
    drawtextcenter!(
        canvas::Canvas, coords::Array{Float64,2}, text::String, color::Color)

Draw characters of an atom and implicit hydrogens between bonds.
"""
function drawtextcenter! end


"""
    drawtextright!(
        canvas::Canvas, coords::Array{Float64,2}, text::String, color::Color)

Draw characters of an atom and implicit hydrogens in a left-to-right direction.
"""
function drawtextright! end


"""
    drawtextannot!(
        canvas::Canvas, coords::Array{Float64,2}, text::String, color::Color,
        bgcolor::Color
    )

Draw a small text next to the atom that shows an atom index or annotation.
"""
function drawtextannot! end


"""
    drawtexthighlight!(canvas::Canvas, coords::Array{Float64,2}, color::Color)

Add a highlight to the atom of interest.
"""
function drawtexthighlight! end


"""
    drawline!(canvas::Canvas, seg::Segment{T<:Point2D}, color::Color)
    drawline!(canvas::Canvas, seg::Segment{T<:Point2D}, ucolor::Color, vcolor::Color)

Draw a solid line.
"""
function drawline! end


"""
    drawdashedline!(canvas::Canvas, seg::Segment{T<:Point2D}, color::Color)
    drawdashedline!(canvas::Canvas, seg::Segment{T<:Point2D}, ucolor::Color, vcolor::Color)

Draw a dashed line.
"""
function drawdashedline! end


"""
    drawwedge!(canvas::Canvas, seg::Segment{T<:Point2D}, color::Color)
    drawwedge!(canvas::Canvas, seg::Segment{T<:Point2D}, ucolor::Color, vcolor::Color)

Draw a wedge.
"""
function drawwedge! end


"""
    drawdashedwedge!(canvas::Canvas, seg::Segment{T<:Point2D}, color::Color)
    drawdashedwedge!(canvas::Canvas, seg::Segment{T<:Point2D}, ucolor::Color, vcolor::Color)

Draw a dashed wedge.
"""
function drawdashedwedge! end


"""
    drawwave!(canvas::Canvas, seg::Segment{T<:Point2D}, color::Color)
    drawwave!(canvas::Canvas, seg::Segment{T<:Point2D}, ucolor::Color, vcolor::Color)

Draw a waved line.
"""
function drawwave! end


"""
    drawlinehighlight!(canvas::Canvas, seg::Segment{T<:Point2D}, color::Color)

Add a highlight to the line.
"""
function drawlinehighlight! end