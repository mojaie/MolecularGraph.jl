#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

module CairoExt


using Colors:
    Colors, RGB, coloralpha

using GeometryBasics:
    GeometryBasics, Point2d, norm

using Graphs:
    Graphs, Edge, src, dst

using MolecularGraph:
    MolecularGraph, Canvas, SimpleMolGraph,
    initcanvas!,
    drawtext!, drawtextannot!, drawtexthighlight!,
    drawline!, drawdashedline!, drawwedge!, drawdashedwedge!,
    drawwave!, drawlinehighlight!,
    is_atom_visible, transformmatrix,
    draw2d!, sethighlight!, drawatomindex!

using Cairo:
    Cairo, CairoContext, CairoMatrix,
    CairoSurface, CairoARGBSurface,
    set_source, get_operator, set_operator,
    get_matrix, set_matrix, set_font_face,
    set_line_width, set_line_cap, set_dash,
    paint, pop_group, restore, save, write_to_png,
    translate, scale,
    move_to, line_to, stroke, close_path,
    arc, rectangle, fill, text, text_extents


mutable struct CairoCanvas <: Canvas
    scaleunit::Float64
    mbwidth::Float64
    wedgewidth::Float64
    wavewidth::Float64
    triminner::Float64
    trimoverlap::Float64
    linehlwidth::Float64
    annotsizef::Float64
    hlsizef::Float64
    paddingXf::Float64
    paddingYf::Float64

    fontweight::String
    fontfamily::String
    fontsize::Float64
    fonttagmap::Dict{Symbol,Tuple{String,String}}
    bgcolor::RGB
    bgopacity::Float64

    cairoscalef::Float64
    surface::Cairo.CairoSurface
    context::Cairo.CairoContext
    coords::Vector{Point2d}

    function CairoCanvas(width::Int, height::Int, bgcolor::RGB, bgopacity::Float64)
        canvas = new()

        # Geometry
        canvas.scaleunit = 30.0
        canvas.mbwidth = 4.5
        canvas.wedgewidth = 4.5
        canvas.wavewidth = 4.5
        canvas.triminner = 3.0
        canvas.trimoverlap = 9.0
        canvas.linehlwidth = 9.0
        canvas.annotsizef = 0.7
        canvas.hlsizef = 0.8
        canvas.paddingXf = 1.0
        canvas.paddingYf = 1.0

        # Appearance
        canvas.fontweight = "Normal"
        canvas.fontfamily = "Sans"
        canvas.fontsize = 11.0
        canvas.fonttagmap = Dict(
            :sub => ("<sub>", "</sub>"),
            :sup => ("<sup>", "</sup>")
        )
        canvas.bgcolor = bgcolor
        canvas.bgopacity = bgopacity

        # Canvas state
        canvas.cairoscalef = 1
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
function MolecularGraph.drawpng(io::IO, mol::SimpleMolGraph, width::Int, height::Int;
        bgcolor="#FFF", bgopacity=1.0,
        atomhighlight=eltype(mol)[], bondhighlight=Edge{eltype(mol)}[], highlightcolor="#FDD835",
        atomindex=false, indexcolor="#000", indexbgcolor="#F0F0FF",
        kwargs...)
    canvas = CairoCanvas(width, height, parse(RGB, bgcolor), bgopacity)
    # draw background
    Cairo.set_source(
        canvas.context, coloralpha(canvas.bgcolor, canvas.bgopacity))
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
    sethighlight!(canvas, nodes_to_show, parse(RGB, highlightcolor))
    sethighlight!(canvas, bondhighlight, parse(RGB, highlightcolor))
    Cairo.set_source(canvas.context, Cairo.pop_group(canvas.context))
    Cairo.paint(canvas.context)
    atomindex && drawatomindex!(canvas, is_atom_visible(mol), parse(RGB, indexcolor), parse(RGB, indexbgcolor))
    Cairo.write_to_png(canvas.surface, io)
    return
end

    
"""
    initcanvas!(canvas::Canvas, coords::AbstractArray{Float64}, boundary::Tuple)

Move and adjust the size of the molecule for drawing.
"""
function MolecularGraph.initcanvas!(
        canvas::CairoCanvas, coords::Vector{Point2d}, boundary::Tuple)
    (top, left, width, height, unit) = boundary
    sf = canvas.scaleunit / unit
    pd = [canvas.paddingXf canvas.paddingYf] * canvas.scaleunit
    conv = p -> (p - Point2d(left, top)) * Point2d(1, -1) * sf + Point2d(pd...)
    canvas.coords = conv.(coords)
    canvasW = canvas.surface.width
    canvasH = canvas.surface.height
    x_scale, y_scale = [canvasW canvasH] ./ (([width height] * sf) .+ (pd * 2))
    min_scale = min(x_scale, y_scale)  # keep aspect ratio
    xoff = x_scale > y_scale ? (1 - y_scale / x_scale) / 2 * canvasW : 0  # x centering
    yoff = y_scale > x_scale ? (1 - x_scale / y_scale) / 2 * canvasH : 0  # y centering
    Cairo.translate(canvas.context, xoff, yoff)
    Cairo.scale(canvas.context, min_scale, min_scale)
    canvas.cairoscalef = min_scale
    return
end


const CAIRO_ATOM_TEXT_HALIGN = Dict(:left => "right", :center => "center", :right => "left")
const CAIRO_ATOM_TEXT_XOFFSET = Dict(:left => 1, :center => 0, :right => -1)


function MolecularGraph.drawtext!(canvas::CairoCanvas, pos, text, color, align)
    halign = CAIRO_ATOM_TEXT_HALIGN[align]
    xoff = CAIRO_ATOM_TEXT_XOFFSET[align]
    fs = round(canvas.fontsize * canvas.cairoscalef, digits=1)
    Cairo.set_font_face(canvas.context, join([canvas.fontfamily, canvas.fontweight, fs], " "))
    Cairo.set_source(canvas.context, color)
    cext = Cairo.text_extents(canvas.context, "C")
    exof = Point2d(0.5, -2.0)  # TODO: extra offset
    xoffset = cext[3] * xoff  # extent of character "C": xb, yb, w, h, xa, ya
    yoffset = cext[4] / 2
    pos_ = pos + Point2d(xoffset, yoffset) + exof
    Cairo.text(
        canvas.context, pos_[1], pos_[2],
        text, halign=halign, valign="center", markup=true
    )
    return
end


function MolecularGraph.drawtextannot!(canvas::CairoCanvas, pos, text, color, bgcolor)
    font_size = round(canvas.fontsize * canvas.annotsizef * canvas.cairoscalef, digits=1)
    bg_size = round(canvas.fontsize * canvas.annotsizef, digits=1)
    Cairo.set_font_face(canvas.context, join([canvas.fontfamily, canvas.fontweight, font_size], " "))
    Cairo.set_source(canvas.context, bgcolor)
    Cairo.arc(canvas.context, pos[1] + bg_size, pos[2] + bg_size, bg_size, 0, 2pi)
    Cairo.fill(canvas.context)
    Cairo.set_source(canvas.context, color)
    Cairo.text(canvas.context, pos[1] + bg_size/2, pos[2] + bg_size/2, text, halign="left", valign="top")
    return
end

function MolecularGraph.drawtexthighlight!(canvas::CairoCanvas, pos, color)
    size = round(Int, canvas.fontsize * canvas.hlsizef)
    Cairo.set_source(canvas.context, color)
    Cairo.arc(canvas.context, pos[1], pos[2], size, 0, 2pi)
    Cairo.fill(canvas.context)
    return
end


function MolecularGraph.drawline!(canvas::CairoCanvas, seg, color; isdashed=false)
    Cairo.set_source(canvas.context, color)
    Cairo.set_line_width(canvas.context, canvas.cairoscalef)
    isdashed && Cairo.set_dash(canvas.context, [10, 10], 0)
    Cairo.move_to(canvas.context, seg[1][1], seg[1][2])
    Cairo.line_to(canvas.context, seg[2][1], seg[2][2])
    Cairo.stroke(canvas.context)
    return
end

function MolecularGraph.drawline!(canvas::CairoCanvas, seg, ucolor, vcolor; isdashed=false)
    ucolor == vcolor && return drawline!(canvas, seg, ucolor, isdashed=isdashed)
    mid = (seg[1] + seg[2]) / 2
    Cairo.set_line_width(canvas.context, canvas.cairoscalef)
    isdashed && Cairo.set_dash(canvas.context, [10, 10], 0)
    Cairo.set_source(canvas.context, ucolor)
    Cairo.move_to(canvas.context, seg[1][1], seg[1][2])
    Cairo.line_to(canvas.context, mid[1], mid[2])
    Cairo.stroke(canvas.context)
    Cairo.set_source(canvas.context, vcolor)
    Cairo.move_to(canvas.context, mid[1], mid[2])
    Cairo.line_to(canvas.context, seg[2][1], seg[2][2])
    Cairo.stroke(canvas.context)
    return
end

MolecularGraph.drawdashedline!(canvas::CairoCanvas, seg, ucolor, vcolor) = drawline!(
    canvas, seg, ucolor, vcolor, isdashed=true)


function cairo_transform!(ctx, transform_mat)
    # Emulate cairo_transform. Be sure to do save-transform.
    m = Cairo.get_matrix(ctx)
    tm = [m.xx m.xy m.x0; m.yx m.yy m.y0; 1 1 1] * transform_mat
    newm = Cairo.CairoMatrix(tm[1, 1], tm[2, 1], tm[1, 2], tm[2, 2], tm[1, 3], tm[2, 3])
    Cairo.set_matrix(ctx, newm)
    return
end

function MolecularGraph.drawwedge!(canvas::CairoCanvas, seg, color)
    """ u ◀︎ v """
    Cairo.save(canvas.context)
    cairo_transform!(
        canvas.context, transformmatrix(seg, 1.0, canvas.wedgewidth))
    Cairo.set_source(canvas.context, color)
    # draw '◀︎', length 1, width 2
    Cairo.move_to(canvas.context, 0, 0)
    Cairo.line_to(canvas.context, 1, 1)
    Cairo.line_to(canvas.context, 1, -1)
    Cairo.close_path(canvas.context)
    Cairo.fill(canvas.context)
    Cairo.restore(canvas.context)
    return
end

function MolecularGraph.drawwedge!(canvas::CairoCanvas, seg, ucolor, vcolor)
    """ u ◀︎ v """
    ucolor == vcolor && return drawwedge!(canvas, seg, ucolor)
    Cairo.save(canvas.context)
    cairo_transform!(
        canvas.context, transformmatrix(seg, 1.0, canvas.wedgewidth))
    Cairo.set_source(canvas.context, ucolor)
    # draw '◀︎', length 1, width 2
    Cairo.move_to(canvas.context, 0, 0)
    Cairo.line_to(canvas.context, 0.5, 0.5)
    Cairo.line_to(canvas.context, 0.5, -0.5)
    Cairo.close_path(canvas.context)
    Cairo.fill(canvas.context)
    Cairo.set_source(canvas.context, vcolor)
    Cairo.move_to(canvas.context, 0.5, 0.5)
    Cairo.line_to(canvas.context, 1, 1)
    Cairo.line_to(canvas.context, 1, -1)
    Cairo.line_to(canvas.context, 0.5, -0.5)
    Cairo.close_path(canvas.context)
    Cairo.fill(canvas.context)
    Cairo.restore(canvas.context)
    return
end


function MolecularGraph.drawdashedwedge!(canvas::CairoCanvas, seg, color)
    """ u ◁ v """
    Cairo.save(canvas.context)
    cairo_transform!(
        canvas.context, transformmatrix(seg, 1 / 7, canvas.wedgewidth / 8))
    Cairo.set_line_width(canvas.context, norm(seg[2] - seg[1]) / 20 * canvas.cairoscalef)  # 16 seems a bit thick
    Cairo.set_source(canvas.context, color)
    # draw '◀︎', length 7, width 16
    for i in 1:8
        Cairo.move_to(canvas.context, i - 1, i)
        Cairo.line_to(canvas.context, i - 1, -i)
        Cairo.stroke(canvas.context)
    end
    Cairo.restore(canvas.context)
    return
end

function MolecularGraph.drawdashedwedge!(canvas::CairoCanvas, seg, ucolor, vcolor)
    """ u ◁ v """
    ucolor == vcolor && return drawdashedwedge!(canvas, seg, ucolor)
    Cairo.save(canvas.context)
    cairo_transform!(
        canvas.context, transformmatrix(seg, 1 / 7, canvas.wedgewidth / 8))
    Cairo.set_line_width(canvas.context, norm(seg[2] - seg[1]) / 20 * canvas.cairoscalef)  # 16 seems a bit thick
    Cairo.set_source(canvas.context, ucolor)
    # draw '◀︎', length 7, width 16
    for i in 1:4
        Cairo.move_to(canvas.context, i - 1, i)
        Cairo.line_to(canvas.context, i - 1, -i)
        Cairo.stroke(canvas.context)
    end
    Cairo.set_source(canvas.context, vcolor)
    for i in 5:8
        Cairo.move_to(canvas.context, i - 1, i)
        Cairo.line_to(canvas.context, i - 1, -i)
        Cairo.stroke(canvas.context)
    end
    Cairo.restore(canvas.context)
    return
end


function MolecularGraph.drawwave!(canvas::CairoCanvas, seg, color)
    Cairo.save(canvas.context)
    cairo_transform!(
        canvas.context, transformmatrix(seg, 1 / 7, canvas.wedgewidth))
    Cairo.set_line_width(canvas.context, canvas.cairoscalef)
    Cairo.set_source(canvas.context, color)
    # draw '~', length 7, width 2
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

function MolecularGraph.drawwave!(canvas::CairoCanvas, seg, ucolor, vcolor)
    ucolor == vcolor && return drawwave!(canvas, seg, ucolor)
    Cairo.save(canvas.context)
    cairo_transform!(
        canvas.context, transformmatrix(seg, 1 / 7, canvas.wedgewidth))
    Cairo.set_line_width(canvas.context, canvas.cairoscalef)
    Cairo.set_source(canvas.context, ucolor)
    # draw '~', length 7, width 2
    Cairo.move_to(canvas.context, 0, 0)
    Cairo.line_to(canvas.context, 0.5, 0)
    Cairo.line_to(canvas.context, 1, 1)
    Cairo.line_to(canvas.context, 2, -1)
    Cairo.line_to(canvas.context, 3, 1)
    Cairo.line_to(canvas.context, 3.5, 0)
    Cairo.stroke(canvas.context)
    Cairo.set_source(canvas.context, vcolor)
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


function MolecularGraph.drawlinehighlight!(canvas::CairoCanvas, seg, color)
    w = round(Int, canvas.linehlwidth * canvas.cairoscalef)
    Cairo.set_source(canvas.context, color)
    Cairo.set_line_cap(canvas.context, Cairo.CAIRO_LINE_CAP_ROUND)
    Cairo.set_line_width(canvas.context, w)
    Cairo.move_to(canvas.context, seg[1][1], seg[1][2])
    Cairo.line_to(canvas.context, seg[2][1], seg[2][2])
    Cairo.stroke(canvas.context)
    return
end


end # module