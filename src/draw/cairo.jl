#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

import Cairo


mutable struct CairoCanvas <: Canvas
    fontweight::String
    fontfamily::String
    fontsize::Float64
    fontscalef::Float64
    bgcolor::Color
    bgopacity::Float64

    mbwidthf::Float64
    wedgewidthf::Float64
    wavewidthf::Float64
    triminnerf::Float64
    trimoverlapf::Float64
    annotsizef::Float64
    scalef::Float64
    paddingX::Float64
    paddingY::Float64

    surface::Cairo.CairoSurface
    context::Cairo.CairoContext
    coords::Matrix{Float64}
    valid::Bool

    function CairoCanvas(width::Int, height::Int, bgcolor::Color, bgopacity::Float64)
        canvas = new()
        canvas.fontweight = "Normal"
        canvas.fontfamily = "Sans"
        canvas.fontsize = 11.0
        canvas.fontscalef = 1
        canvas.bgcolor = bgcolor
        canvas.bgopacity = bgopacity

        canvas.mbwidthf = 0.15
        canvas.wedgewidthf = 0.3
        canvas.wavewidthf = 0.2
        canvas.triminnerf = 0.2
        canvas.trimoverlapf = 0.3
        canvas.annotsizef = 0.7
        canvas.scalef = 30.0 # suitable for fontsize=11 and default line width
        canvas.paddingX = 30.0
        canvas.paddingY = 30.0

        canvas.surface = Cairo.CairoARGBSurface(width, height)
        canvas.context = Cairo.CairoContext(canvas.surface)

        return canvas
    end
end



"""
    drawpng(io::IO, mol::SimpleMolGraph, width::Int, height::Int; kwargs...)

Generate molecular structure image as a PNG format.

`width` and `height` specifies the size of the image in px.
"""
function drawpng(io::IO, mol::SimpleMolGraph, width::Int, height::Int;
        bgcolor=Color(255.0, 255.0, 255.0), bgopacity=0.0,
        atomhighlight=eltype(mol)[], bondhighlight=Edge{eltype(mol)}[], highlightcolor=Color(253, 216, 53),
        atomindex=false, indexcolor=Color(0, 0, 0), indexbgcolor=Color(240, 240, 255),
        kwargs...)
    canvas = CairoCanvas(width, height, bgcolor, bgopacity)
    # draw background
    Cairo.set_source_rgba(
        canvas.context, canvas.bgcolor.r / 255, canvas.bgcolor.g / 255, canvas.bgcolor.b / 255, canvas.bgopacity)
    Cairo.rectangle(canvas.context, 0.0, 0.0, canvas.surface.width, canvas.surface.height)  
    Cairo.fill(canvas.context)
    # draw molecule
    draw2d!(canvas, mol; kwargs...)
    # highlight atoms if is_atom_visible=true or no incident edges
    # setdiff(Int[], []) -> Any[], setdiff(Int[], Int[]) -> Int[]  ???
    enodes = Set{eltype(mol)}(vcat([[src(e), dst(e)] for e in bondhighlight]...))
	nodes_to_show = collect(setdiff(
        atomhighlight, setdiff(enodes, findall(is_atom_visible(mol)))))
    # set highlights behind the elements.
    op = Cairo.get_operator(canvas.context)  # Cairo.OPERATOR_OVER
    Cairo.set_operator(canvas.context, Cairo.OPERATOR_MULTIPLY)  
    Cairo.push_group(canvas.context)
    Cairo.set_operator(canvas.context, op)  # avoid blending overlapped highlights
    sethighlight!(canvas, nodes_to_show, highlightcolor)
    sethighlight!(canvas, bondhighlight, highlightcolor)
    Cairo.set_source(canvas.context, Cairo.pop_group(canvas.context))
    Cairo.paint(canvas.context)
    atomindex && drawatomindex!(canvas, is_atom_visible(mol), indexcolor, indexbgcolor)
    Cairo.write_to_png(canvas.surface, io)
    return
end

    
"""
    initcanvas!(canvas::Canvas, coords::AbstractArray{Float64}, boundary::Tuple)

Move and adjust the size of the molecule for drawing.
"""
function initcanvas!(
        canvas::CairoCanvas, coords::AbstractArray{Float64}, boundary::Tuple)
    (top, left, width, height, unit) = boundary
    sf = canvas.scalef / unit
    pd = [canvas.paddingX canvas.paddingY]
    canvas.coords = (coords .- [left top]) .* [1 -1] .* sf .+ pd
    viewboxW = width * sf + canvas.paddingX * 2
    viewboxH = height * sf + canvas.paddingY * 2
    x_scale = canvas.surface.width / viewboxW
    y_scale = canvas.surface.height / viewboxH
    min_scale = min(x_scale, y_scale)  # keep aspect ratio
    xoff = x_scale > y_scale ? (1 - y_scale / x_scale) / 2 * canvas.surface.width : 0  # x centering
    yoff = y_scale > x_scale ? (1 - x_scale / y_scale) / 2 * canvas.surface.height : 0  # y centering
    Cairo.translate(canvas.context, xoff, yoff)
    Cairo.scale(canvas.context, min_scale, min_scale)
    canvas.fontscalef = min_scale
    canvas.valid = true
    return
end


atommarkupleft(canvas::CairoCanvas, symbol, charge, implicith) = atomhtml(
    symbol, charge, implicith, :left)

atommarkupright(canvas::CairoCanvas, symbol, charge, implicith) = atomhtml(
    symbol, charge, implicith, :right)


function drawtextcairo!(canvas::CairoCanvas, pos, text, color, fxoff, halign)
    fs = round(canvas.fontsize * canvas.fontscalef, digits=1)
    Cairo.set_font_face(canvas.context, join([canvas.fontfamily, canvas.fontweight, fs], " "))
    Cairo.set_source_rgba(canvas.context, color.r / 255, color.g / 255, color.b / 255, 1)
    cext = Cairo.text_extents(canvas.context, "C")
    xoffset = fxoff(cext[3])  # extent of character "C": xb, yb, w, h, xa, ya
    yoffset = cext[4] / 2
    Cairo.text(
        canvas.context, pos.x + xoffset, pos.y + yoffset,
        text, halign=halign, valign="center", markup=true
    )
    return
end

drawtextleft!(canvas::CairoCanvas, pos, text, color) = drawtextcairo!(
    canvas, pos, text, color, x -> x, "right")

drawtextcenter!(canvas::CairoCanvas, pos, text, color) = drawtextcairo!(
    canvas, pos, text, color, x -> 0, "center")

drawtextright!(canvas::CairoCanvas, pos, text, color) = drawtextcairo!(
    canvas, pos, text, color, x -> -x, "left")


function drawtextannot!(canvas::CairoCanvas, pos, text, color, bgcolor)
    size = round(canvas.fontsize * canvas.fontscalef * canvas.annotsizef, digits=1)
    Cairo.set_font_face(canvas.context, join([canvas.fontfamily, canvas.fontweight, size], " "))
    Cairo.set_source_rgba(canvas.context, bgcolor.r / 255, bgcolor.g / 255, bgcolor.b / 255, 1)
    Cairo.arc(canvas.context, pos.x + size, pos.y + size, size, 0, 2pi)
    Cairo.fill(canvas.context)
    Cairo.set_source_rgba(canvas.context, color.r / 255, color.g / 255, color.b / 255, 1)
    Cairo.text(canvas.context, pos.x, pos.y, text, halign="left", valign="top")
    return
end

function drawtexthighlight!(canvas::CairoCanvas, pos, color)
    size = round(Int, canvas.fontsize * canvas.fontscalef * 1.2)
    Cairo.set_source_rgba(canvas.context, color.r / 255, color.g / 255, color.b / 255, 1)
    Cairo.arc(canvas.context, pos.x, pos.y, size, 0, 2pi)
    Cairo.fill(canvas.context)
    return
end


function drawline!(canvas::CairoCanvas, seg, color; isdashed=false)
    Cairo.set_source_rgba(canvas.context, color.r / 255, color.g / 255, color.b / 255, 1)
    Cairo.set_line_width(canvas.context, canvas.fontscalef)
    isdashed && Cairo.set_dash(canvas.context, [10, 10], 0)
    Cairo.move_to(canvas.context, seg.u.x, seg.u.y)
    Cairo.line_to(canvas.context, seg.v.x, seg.v.y)
    Cairo.stroke(canvas.context)
    return
end

function drawline!(canvas::CairoCanvas, seg, ucolor, vcolor; isdashed=false)
    ucolor == vcolor && return drawline!(canvas, seg, ucolor, isdashed=isdashed)
    mid = midpoint(seg)
    Cairo.set_line_width(canvas.context, canvas.fontscalef)
    isdashed && Cairo.set_dash(canvas.context, [10, 10], 0)
    Cairo.set_source_rgba(canvas.context, ucolor.r / 255, ucolor.g / 255, ucolor.b / 255, 1)
    Cairo.move_to(canvas.context, seg.u.x, seg.u.y)
    Cairo.line_to(canvas.context, mid.x, mid.y)
    Cairo.stroke(canvas.context)
    Cairo.set_source_rgba(canvas.context, vcolor.r / 255, vcolor.g / 255, vcolor.b / 255, 1)
    Cairo.move_to(canvas.context, mid.x, mid.y)
    Cairo.line_to(canvas.context, seg.v.x, seg.v.y)
    Cairo.stroke(canvas.context)
    return
end

drawdashedline!(canvas::CairoCanvas, seg, ucolor, vcolor) = drawline!(
    canvas, seg, ucolor, vcolor, isdashed=true)


function cairo_transform!(ctx, scalef, rotatef, translf)
    # Emulate cairo_transform. Be sure to do save-transform.
    m = Cairo.get_matrix(ctx)
    tm = [m.xx m.xy m.x0; m.yx m.yy m.y0; 1 1 1] * transformmatrix(scalef, rotatef, translf)
    newm = Cairo.CairoMatrix(tm[1, 1], tm[2, 1], tm[1, 2], tm[2, 2], tm[1, 3], tm[2, 3])
    Cairo.set_matrix(ctx, newm)
    return
end

function drawwedge!(canvas::CairoCanvas, seg, color)
    """ u ◀︎ v """
    d = distance(seg)
    scalef = Point2D(d, d / 2 * canvas.wedgewidthf)
    rotatef = unitvector(seg)
    translf = seg.u
    Cairo.save(canvas.context)
    cairo_transform!(canvas.context, scalef, rotatef, translf)
    Cairo.set_source_rgba(canvas.context, color.r / 255, color.g / 255, color.b / 255, 1)
    Cairo.move_to(canvas.context, 0, 0)
    Cairo.line_to(canvas.context, 1, 1)
    Cairo.line_to(canvas.context, 1, -1)
    Cairo.close_path(canvas.context)
    Cairo.fill(canvas.context)
    Cairo.restore(canvas.context)
    return
end

function drawwedge!(canvas::CairoCanvas, seg, ucolor, vcolor)
    """ u ◀︎ v """
    ucolor == vcolor && return drawwedge!(canvas, seg, ucolor)
    d = distance(seg)
    scalef = Point2D(d, d / 2 * canvas.wedgewidthf)
    rotatef = unitvector(seg)
    translf = seg.u
    Cairo.save(canvas.context)
    cairo_transform!(canvas.context, scalef, rotatef, translf)
    Cairo.set_source_rgba(canvas.context, ucolor.r / 255, ucolor.g / 255, ucolor.b / 255, 1)
    Cairo.move_to(canvas.context, 0, 0)
    Cairo.line_to(canvas.context, 0.5, 0.5)
    Cairo.line_to(canvas.context, 0.5, -0.5)
    Cairo.close_path(canvas.context)
    Cairo.fill(canvas.context)
    Cairo.set_source_rgba(canvas.context, vcolor.r / 255, vcolor.g / 255, vcolor.b / 255, 1)
    Cairo.move_to(canvas.context, 0.5, 0.5)
    Cairo.line_to(canvas.context, 1, 1)
    Cairo.line_to(canvas.context, 1, -1)
    Cairo.line_to(canvas.context, 0.5, -0.5)
    Cairo.close_path(canvas.context)
    Cairo.fill(canvas.context)
    Cairo.restore(canvas.context)
    return
end


function drawdashedwedge!(canvas::CairoCanvas, seg, color)
    """ u ◁ v """
    d = distance(seg)
    scalef = Point2D(d / 7, d / 14 * canvas.wedgewidthf)
    rotatef = unitvector(seg)
    translf = seg.u
    Cairo.save(canvas.context)
    cairo_transform!(canvas.context, scalef, rotatef, translf)
    Cairo.set_line_width(canvas.context, d / 20 * canvas.fontscalef)  # 16 seems a bit thick
    Cairo.set_source_rgba(canvas.context, color.r / 255, color.g / 255, color.b / 255, 1)
    for i in 1:8
        Cairo.move_to(canvas.context, i - 1, i)
        Cairo.line_to(canvas.context, i - 1, -i)
        Cairo.stroke(canvas.context)
    end
    Cairo.restore(canvas.context)
    return
end

function drawdashedwedge!(canvas::CairoCanvas, seg, ucolor, vcolor)
    """ u ◁ v """
    ucolor == vcolor && return drawdashedwedge!(canvas, seg, ucolor)
    d = distance(seg)
    scalef = Point2D(d / 7, d / 14 * canvas.wedgewidthf)
    rotatef = unitvector(seg)
    translf = seg.u
    Cairo.save(canvas.context)
    cairo_transform!(canvas.context, scalef, rotatef, translf)
    Cairo.set_line_width(canvas.context, d / 20 * canvas.fontscalef)  # 16 seems a bit thick
    Cairo.set_source_rgba(canvas.context, ucolor.r / 255, ucolor.g / 255, ucolor.b / 255, 1)
    for i in 1:4
        Cairo.move_to(canvas.context, i - 1, i)
        Cairo.line_to(canvas.context, i - 1, -i)
        Cairo.stroke(canvas.context)
    end
    Cairo.set_source_rgba(canvas.context, vcolor.r / 255, vcolor.g / 255, vcolor.b / 255, 1)
    for i in 5:8
        Cairo.move_to(canvas.context, i - 1, i)
        Cairo.line_to(canvas.context, i - 1, -i)
        Cairo.stroke(canvas.context)
    end
    Cairo.restore(canvas.context)
    return
end


function drawwave!(canvas::CairoCanvas, seg, color)
    d = distance(seg)
    scalef = Point2D(d / 7, d / 2 * canvas.wedgewidthf)
    rotatef = unitvector(seg)
    translf = seg.u
    Cairo.save(canvas.context)
    cairo_transform!(canvas.context, scalef, rotatef, translf)
    Cairo.set_line_width(canvas.context, canvas.fontscalef)
    Cairo.set_source_rgba(canvas.context, color.r / 255, color.g / 255, color.b / 255, 1)
    Cairo.move_to(canvas.context, 0, 0)
    Cairo.line_to(canvas.context, 0.5, 0)
    Cairo.line_to(canvas.context, 1, 1)
    Cairo.line_to(canvas.context, 2, -1)
    Cairo.line_to(canvas.context, 3, 1)
    Cairo.line_to(canvas.context, 4, -1)
    Cairo.line_to(canvas.context, 5, 1)
    Cairo.line_to(canvas.context, 6, -1)
    Cairo.line_to(canvas.context, 6.5, 0)
    Cairo.line_to(canvas.context, 7, 0)
    Cairo.stroke(canvas.context)
    Cairo.restore(canvas.context)
    return
end

function drawwave!(canvas::CairoCanvas, seg, ucolor, vcolor)
    ucolor == vcolor && return drawwave!(canvas, seg, ucolor)
    d = distance(seg)
    scalef = Point2D(d / 7, d / 2 * canvas.wedgewidthf)
    rotatef = unitvector(seg)
    translf = seg.u
    Cairo.save(canvas.context)
    cairo_transform!(canvas.context, scalef, rotatef, translf)
    Cairo.set_line_width(canvas.context, canvas.fontscalef)
    Cairo.set_source_rgba(canvas.context, ucolor.r / 255, ucolor.g / 255, ucolor.b / 255, 1)
    Cairo.move_to(canvas.context, 0, 0)
    Cairo.line_to(canvas.context, 0.5, 0)
    Cairo.line_to(canvas.context, 1, 1)
    Cairo.line_to(canvas.context, 2, -1)
    Cairo.line_to(canvas.context, 3, 1)
    Cairo.line_to(canvas.context, 3.5, 0)
    Cairo.stroke(canvas.context)
    Cairo.set_source_rgba(canvas.context, vcolor.r / 255, vcolor.g / 255, vcolor.b / 255, 1)
    Cairo.move_to(canvas.context, 3.5, 0)
    Cairo.line_to(canvas.context, 4, -1)
    Cairo.line_to(canvas.context, 5, 1)
    Cairo.line_to(canvas.context, 6, -1)
    Cairo.line_to(canvas.context, 6.5, 0)
    Cairo.line_to(canvas.context, 7, 0)
    Cairo.stroke(canvas.context)
    Cairo.restore(canvas.context)
    return
end


function drawlinehighlight!(canvas::CairoCanvas, seg, color)
    w = round(Int, 10 * canvas.fontscalef)
    Cairo.set_source_rgba(canvas.context, color.r / 255, color.g / 255, color.b / 255, 1)
    Cairo.set_line_cap(canvas.context, Cairo.CAIRO_LINE_CAP_ROUND)
    Cairo.set_line_width(canvas.context, w)
    Cairo.move_to(canvas.context, seg.u.x, seg.u.y)
    Cairo.line_to(canvas.context, seg.v.x, seg.v.y)
    Cairo.stroke(canvas.context)
    return
end