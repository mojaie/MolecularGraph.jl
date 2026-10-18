#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#


abstract type Canvas end


# Custom 2D drawing systems using `draw2d!` should implement following methods (see ./svg.jl and ./cairo.jl)

"""
    initcanvas!(canvas::Canvas, coords::Vector{Point2d}, width::Float64, height::Float64)

Convert the coordinate system of SDFile to that of the graphics library and set the drawing area.
"""
function initcanvas! end


"""
    drawtext!(
        canvas::Canvas, pos::Point2d, text::String, color::Color, align::Symbol)

Draw characters of an atom and implicit hydrogens in a right-to-left direction.
"""
function drawtext! end


"""
    drawtextannot!(
        canvas::Canvas, pos::Point2d, text::String, color::Color,
        bgcolor::Color
    )

Draw a small text next to the atom that shows an atom index or annotation.
"""
function drawtextannot! end


"""
    drawtexthighlight!(canvas::Canvas, pos::Point2d, color::Color)

Add a highlight to the atom of interest.
"""
function drawtexthighlight! end


"""
    drawline!(canvas::Canvas, seg::Line{2,Float64}, color::Color)
    drawline!(canvas::Canvas, seg::Line{2,Float64}, ucolor::Color, vcolor::Color)

Draw a solid line.
"""
function drawline! end


"""
    drawdashedline!(canvas::Canvas, seg::Line{2,Float64}, color::Color)
    drawdashedline!(canvas::Canvas, seg::Line{2,Float64}, ucolor::Color, vcolor::Color)

Draw a dashed line.
"""
function drawdashedline! end


"""
    drawwedge!(canvas::Canvas, seg::Line{2,Float64}, color::Color)
    drawwedge!(canvas::Canvas, seg::Line{2,Float64}, ucolor::Color, vcolor::Color)

Draw a wedge.
"""
function drawwedge! end


"""
    drawdashedwedge!(canvas::Canvas, seg::Line{2,Float64}, color::Color)
    drawdashedwedge!(canvas::Canvas, seg::Line{2,Float64}, ucolor::Color, vcolor::Color)

Draw a dashed wedge.
"""
function drawdashedwedge! end


"""
    drawwave!(canvas::Canvas, seg::Line{2,Float64}, color::Color)
    drawwave!(canvas::Canvas, seg::Line{2,Float64}, ucolor::Color, vcolor::Color)

Draw a waved line.
"""
function drawwave! end


"""
    drawlinehighlight!(canvas::Canvas, seg::Line{2,Float64}, color::Color)

Add a highlight to the line.
"""
function drawlinehighlight! end


"""
    drawpng(io::IO, mol::SimpleMolGraph, width::Int, height::Int; kwargs...)

Generate molecular structure image as a PNG format.

`width` and `height` specifies the size of the image in px.
"""
function drawpng end


"""
    spacefilling(mol::MolGraph; radii="van der Waals")

Represent `mol` as a space-filling (Calotte) model in three dimensions. `mol` should have 3d atom
positions in Angstroms. Requires loading a Makie backend (e.g. GLMakie or CairoMakie).
"""
function spacefilling end

"""
    spacefilling!(ax, mol::MolGraph; radii="van der Waals")

In-place version of [`spacefilling`](@ref) that adds the model to an existing Makie axis or scene `ax`.
"""
function spacefilling! end

"""
    ballstick(mol::MolGraph; radii=0.4, bonddiameter=0.1)

Represent `mol` as a ball-and-stick model in three dimensions. `mol` should have 3d atom positions
in Angstroms. Requires loading a Makie backend (e.g. GLMakie or CairoMakie).
"""
function ballstick end

"""
    ballstick!(ax, mol::MolGraph; radii=0.4, bonddiameter=0.1)

In-place version of [`ballstick`](@ref) that adds the model to an existing Makie axis or scene `ax`.
"""
function ballstick! end

"""
    stick(mol::MolGraph; size=0.33)

Represent `mol` as a stick model in three dimensions. `mol` should have 3d atom positions in Angstroms.
Requires loading a Makie backend (e.g. GLMakie or CairoMakie).
"""
function stick end

"""
    stick!(ax, mol::MolGraph; size=0.33)

In-place version of [`stick`](@ref) that adds the model to an existing Makie axis or scene `ax`.
"""
function stick! end

"""
    wire(mol::MolGraph; size=0.1)

Represent `mol` as a wire-frame model in three dimensions. `mol` should have 3d atom positions in Angstroms.
Requires loading a Makie backend (e.g. GLMakie or CairoMakie).
"""
function wire end

"""
    wire!(ax, mol::MolGraph; size=0.1)

In-place version of [`wire`](@ref) that adds the model to an existing Makie axis or scene `ax`.
"""
function wire! end