#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    SvgCanvas,
    tosvg,
    drawsvg,
    initcanvas!,
    savesvg

using Printf


mutable struct SvgCanvas <: Canvas
    fontweight::String
    fontfamily::String
    fontsize::Float64
    bgcolor::Color
    opacity::Float64

    mbwidthf::Float64
    wedgewidthf::Float64
    wavewidthf::Float64
    triminnerf::Float64
    trimoverlapf::Float64
    annotsizef::Float64
    scalef::Float64
    paddingX::Float64
    paddingY::Float64

    viewboxW::Float64
    viewboxH::Float64
    
    bgelements::Vector{String}
    elements::Vector{String}
    coords::Matrix{Float64}
    valid::Bool

    function SvgCanvas()
        canvas = new()
        canvas.fontweight = "normal"
        canvas.fontfamily = "Helvetica"
        canvas.fontsize = 14.0
        canvas.bgcolor = Color(255.0, 255.0, 255.0)
        canvas.opacity = 0.0

        canvas.mbwidthf = 0.15
        canvas.wedgewidthf = 0.3
        canvas.wavewidthf = 0.2
        canvas.triminnerf = 0.2
        canvas.trimoverlapf = 0.3
        canvas.annotsizef = 0.7
        canvas.scalef = 30.0 # suitable for fontsize=14 and default line width
        canvas.paddingX = 30.0
        canvas.paddingY = 30.0

        canvas.bgelements = []
        canvas.elements = []
        canvas.valid = false
        return canvas
    end
end


svgcolor(c::Color) = @sprintf "rgb(%d, %d, %d)" c.r c.g c.b
svgcoords(p::Point2D) = @sprintf "x=\"%.2f\" y=\"%.2f\"" p.x p.y
svgcirclecoords(p::Point2D) = @sprintf "cx=\"%.2f\" cy=\"%.2f\"" p.x p.y
svgcoords(s::Segment
    ) = @sprintf "x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\"" s.u.x s.u.y s.v.x s.v.y
svgtransform(tf::Array{Float64,2}
    ) = @sprintf "%.2f %.2f %.2f %.2f %.2f %.2f" tf[1] tf[2] tf[4] tf[5] tf[7] tf[8]


function tosvg(canvas::SvgCanvas, width, height)
    vbWf = @sprintf "%.2f" canvas.viewboxW
    vbHf = @sprintf "%.2f" canvas.viewboxH
    bgc = svgcolor(canvas.bgcolor)
    header = """<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:xlink="http://www.w3.org/1999/xlink"
     version="1.2" baseProfile="tiny"
     text-rendering="geometricPrecision"
     preserveAspectRatio="xMidYMid meet"
     font-weight="$(canvas.fontweight)"
     font-family="$(canvas.fontfamily)"
     width="$(round(Int, width))" height="$(round(Int, height))"
     viewBox="0 0 $(vbWf) $(vbHf)">
    """
    bg = """<rect x="0" y="0" width="$(vbWf)" height="$(vbHf)"
     fill="$(bgc)" opacity="$(canvas.opacity)"/>
    """
    endsvg = "</svg>\n"
    return join([header, bg, canvas.bgelements..., canvas.elements..., endsvg], "")
end


"""
    drawsvg(mol::GraphMol, width::Int, height::Int)

Generate molecular structure image as a SVG format string.

`width` and `height` specifies the size of the image (width and height
attribute of svg tag).
"""
function drawsvg(mol::UndirectedGraph, width, height; kwargs...)
    canvas = SvgCanvas()
    draw2d!(canvas, mol; kwargs...)
    if haskey(kwargs, :highlight)
        sethighlight!(canvas, kwargs[:highlight])
    end
    if haskey(kwargs, :atomindex) && kwargs[:atomindex]
        drawatomindex!(canvas, mol)
    end
    return tosvg(canvas, width, height)
end


"""
    boundary(mol::GraphMol, coords::AbstractArray{Float64}
        ) -> (top, left, width, height, unit)

Get boundaries and an appropriate bond length unit for the molecule drawing
canvas.
"""
function boundary(mol::GraphMol, coords::AbstractArray{Float64})
    (left, right) = extrema(x_components(coords))
    (bottom, top) = extrema(y_components(coords))
    width = right - left
    height = top - bottom
    dists = []
    # Size unit
    for (u, v) in edgesiter(mol)
        d = Geometry.distance(Point2D(coords, u), Point2D(coords, v))
        if d > 0.0001  # Remove overlapped
            push!(dists, d)
        end
    end
    if isempty(dists)
        long = max(width, height)
        unit = long > 0.0001 ? long / sqrt(nodecount(mol)) : 1
    else
        unit = median(dists) # Median bond length
    end
    return (top, left, width, height, unit)
end

boundary(view::SubgraphView, coords::AbstractArray{Float64}
    ) = boundary(view.graph, coords)


"""
    initcanvas!(canvas::Canvas, coords::AbstractArray{Float64}, boundary::Tuple)

Move and adjust the size of the molecule for drawing.
"""
function initcanvas!(
        canvas::SvgCanvas, coords::AbstractArray{Float64}, boundary::Tuple)
    isempty(coords) && return
    (top, left, width, height, unit) = boundary
    sf = canvas.scalef / unit
    pd = [canvas.paddingX canvas.paddingY]
    canvas.coords = (coords .- [left top]) .* [1 -1] .* sf .+ pd
    canvas.viewboxW = width * sf + canvas.paddingX * 2
    canvas.viewboxH = height * sf + canvas.paddingY * 2
    canvas.valid = true
    return
end



function trimbond(canvas::SvgCanvas, seg, uvis, vvis)
    if uvis && vvis
        return trim_uv(seg, canvas.trimoverlapf * 2)
    elseif uvis
        return trim_u(seg, canvas.trimoverlapf)
    elseif vvis
        return trim_v(seg, canvas.trimoverlapf)
    end
    return seg
end


function singlebond!(canvas::SvgCanvas, seg, ucolor, vcolor, uvis, vvis)
    seg_ = trimbond(canvas, seg, uvis, vvis)
    drawline!(canvas, seg_, ucolor, vcolor)
    return
end


function wedged!(canvas::SvgCanvas, seg, ucolor, vcolor, uvis, vvis)
    seg_ = trimbond(canvas, seg, uvis, vvis)
    drawwedge!(canvas, seg_, ucolor, vcolor)
    return
end


function dashedwedged!(canvas::SvgCanvas, seg, ucolor, vcolor, uvis, vvis)
    seg_ = trimbond(canvas, seg, uvis, vvis)
    drawdashedwedge!(canvas, seg_, ucolor, vcolor)
    return
end


function wavesingle!(canvas::SvgCanvas, seg, ucolor, vcolor, uvis, vvis)
    seg_ = trimbond(canvas, seg, uvis, vvis)
    drawwave!(canvas, seg_, ucolor, vcolor)
    return
end


function doublebond!(canvas::SvgCanvas, seg, ucolor, vcolor, uvis, vvis)
    dist = canvas.scalef * canvas.mbwidthf / 2
    seg_ = trimbond(canvas, seg, uvis, vvis)
    seg1 = translate(seg_, pi / 2, dist)
    seg2 = translate(seg_, -pi / 2, dist)
    drawline!(canvas, seg1, ucolor, vcolor)
    drawline!(canvas, seg2, ucolor, vcolor)
    return
end


function crossdouble!(canvas::SvgCanvas, seg, ucolor, vcolor, uvis, vvis)
    dist = canvas.scalef * canvas.mbwidthf / 2
    seg_ = trimbond(canvas, seg, uvis, vvis)
    cw = translate(seg_, pi / 2, dist)
    ccw = translate(seg_, -pi / 2, dist)
    drawline!(canvas, Segment(cw.u, ccw.v), ucolor, vcolor)
    drawline!(canvas, Segment(ccw.u, cw.v), ucolor, vcolor)
    return
end


function ringdouble!(canvas::SvgCanvas, seg, ucolor, vcolor, uvis, vvis, direction)
    seg_ = trimbond(canvas, seg, uvis, vvis)
    dist = canvas.scalef * canvas.mbwidthf
    segin = trim_uv(translate(seg_, direction, dist), canvas.triminnerf)
    drawline!(canvas, seg_, ucolor, vcolor)
    drawline!(canvas, segin, ucolor, vcolor)
    return
end

# NOTE: the direction is reversed due to x-axis reflection
clockwisedouble!(canvas::SvgCanvas, seg, ucolor, vcolor, uvis, vvis
    ) = ringdouble!(canvas, seg, ucolor, vcolor, uvis, vvis, pi / 2)
counterdouble!(canvas::SvgCanvas, seg, ucolor, vcolor, uvis, vvis
    ) = ringdouble!(canvas, seg, ucolor, vcolor, uvis, vvis, -pi / 2)


function triplebond!(canvas::SvgCanvas, seg, ucolor, vcolor, uvis, vvis)
    dist = canvas.scalef * canvas.mbwidthf
    seg_ = trimbond(canvas, seg, uvis, vvis)
    seg1 = translate(seg_, pi / 2, dist)
    seg2 = translate(seg_, -pi / 2, dist)
    drawline!(canvas, seg_, ucolor, vcolor)
    drawline!(canvas, seg1, ucolor, vcolor)
    drawline!(canvas, seg2, ucolor, vcolor)
    return
end


const BOND_DRAWER = Dict(
    1 => Dict(
        0 => singlebond!,
        1 => wedged!,
        6 => dashedwedged!,
        4 => wavesingle!
    ),
    2 => Dict(
        0 => clockwisedouble!,
        1 => counterdouble!,
        2 => doublebond!,
        3 => crossdouble!
    ),
    3 => Dict(
        0 => triplebond!
    )
)

setbond!(
    canvas::SvgCanvas, order, notation, seg, ucolor, vcolor, uvis, vvis
) = BOND_DRAWER[order][notation](canvas, seg, ucolor, vcolor, uvis, vvis)


function setbondhighlight!(canvas::SvgCanvas, seg, color)
    cds = svgcoords(seg)
    c = svgcolor(color)
    elem = """<line $(cds) stroke="$(c)" stroke-width="10" stroke-linecap="round"/>
    """
    push!(canvas.bgelements, elem)
    return
end



function atomsymbol!(canvas::SvgCanvas, pos, atomsymbol, color, implicith, charge;
                     direction=:right, anchor=""" text-anchor="end" """,
                     xoffset=0)
    xy = svgcoords(pos + (xoffset, canvas.fontsize/2))
    small = round(Int, canvas.fontsize * 0.7)
    content = atommarkup(
        atomsymbol, charge, implicith, direction,
        """<tspan baseline-shift="-25%" font-size="$(small)">""", "</tspan>",
        """<tspan baseline-shift="50%" font-size="$(small)">""", "</tspan>"
    )
    c = svgcolor(color)
    elem = """<text $(xy) font-size="$(canvas.fontsize)"
     fill="$(c)"$(anchor)>$(content)</text>
    """
    push!(canvas.elements, elem)
    return
end

setatomright!(canvas::SvgCanvas, pos, sym, color, hcnt, charge) = atomsymbol!(
    canvas, pos, sym, color, hcnt, charge,
    direction=:left, xoffset=canvas.fontsize/2
)

setatomcenter!(canvas::SvgCanvas, pos, sym, color, hcnt, charge) = atomsymbol!(
    canvas, pos, sym, color, hcnt, charge, anchor=""" text-anchor="middle" """)

setatomleft!(canvas::SvgCanvas, pos, sym, color, hcnt, charge) = atomsymbol!(
    canvas, pos, sym, color, hcnt, charge,
    anchor=" ", xoffset=canvas.fontsize/-2
)


function setatomnote!(canvas::SvgCanvas, pos, text, color, bgcolor, opacity)
    size = round(Int, canvas.fontsize * canvas.annotsizef)
    bxy = svgcoords(pos)
    txy = svgcoords(pos + (0, size))
    c = svgcolor(color)
    bc = svgcolor(bgcolor)
    elem = """<g>
     <rect $(bxy) width="$(size)" height="$(size)" rx="$(size/2)" ry="$(size/2)" fill="$(bc)" opacity="$(opacity)" />
     <text $(txy) font-size="$(size)" fill="$(c)">$(text) </text>
    </g>
    """
    push!(canvas.elements, elem)
    return
end

function setatomhighlight!(canvas::SvgCanvas, pos, color)
    size = round(Int, canvas.fontsize * 1.2)
    xy = svgcoords(pos - (size / 2, size / 2))
    c = svgcolor(color)
    elem = """
     <rect $(xy) width="$(size)" height="$(size)" rx="$(size/2)" ry="$(size/2)" fill="$(c)" />
    """
    push!(canvas.bgelements, elem)
    return
end


function drawline!(canvas::SvgCanvas, seg, color; isdashed=false)
    option = isdashed ? """ stroke-dasharray="10,10" """ : " "
    coords = svgcoords(seg)
    c = svgcolor(color)
    elem = """<line $(coords) stroke="$(c)"$(option)/>"""
    push!(canvas.elements, elem)
    return
end

function drawline!(canvas::SvgCanvas, seg, ucolor, vcolor; isdashed=false)
    ucolor == vcolor && return drawline!(canvas, seg, ucolor, isdashed=isdashed)
    option = isdashed ? """ stroke-dasharray="10,10" """ : " "
    mid = midpoint(seg)
    coords1 = svgcoords(Segment(seg.u, mid))
    coords2 = svgcoords(Segment(mid, seg.v))
    uc = svgcolor(ucolor)
    vc = svgcolor(vcolor)
    elem = """<line $(coords1) stroke="$(uc)"$(option)/>
    <line $(coords2) stroke="$(vc)"$(option)/>
    """
    push!(canvas.elements, elem)
    return
end

drawdashedline!(canvas::SvgCanvas, seg, ucolor, vcolor) = drawline!(
    canvas, seg, ucolor, vcolor, isdashed=true)


function drawwedge!(canvas::SvgCanvas, seg, color)
    """ u ◀︎ v """
    d = Geometry.distance(seg)
    scalef = Point2D(d, d / 2 * canvas.wedgewidthf)
    rotatef = unitvector(seg)
    translf = seg.u
    svgtf = svgtransform(transformmatrix(scalef, rotatef, translf))
    c = svgcolor(color)
    elem = """<polygon points="0,0 1,1 1,-1" fill="$(c)" transform="matrix($(svgtf))"/>
    """
    push!(canvas.elements, elem)
    return
end

function drawwedge!(canvas::SvgCanvas, seg, ucolor, vcolor)
    """ u ◀︎ v """
    ucolor == vcolor && return drawwedge!(canvas, seg, ucolor)
    d = Geometry.distance(seg)
    scalef = Point2D(d, d / 2 * canvas.wedgewidthf)
    rotatef = unitvector(seg)
    translf = seg.u
    svgtf = svgtransform(transformmatrix(scalef, rotatef, translf))
    uc = svgcolor(ucolor)
    vc = svgcolor(vcolor)
    elem = """<g stroke-width="0.3" transform="matrix($(svgtf))">
     <polygon points="0,0 0.5,-0.5 0.5,0.5" fill="$(uc)"/>
     <polygon points="0.5,-0.5 0.5,0.5 1,1 1,-1" fill="$(vc)"/>
     </g>
    """
    push!(canvas.elements, elem)
    return
end


function drawdashedwedge!(canvas::SvgCanvas, seg, color)
    """ u ◁ v """
    d = Geometry.distance(seg)
    scalef = Point2D(d / 7, d / 14 * canvas.wedgewidthf)
    rotatef = unitvector(seg)
    translf = seg.u
    svgtf = svgtransform(transformmatrix(scalef, rotatef, translf))
    c = svgcolor(color)
    elem = """<g stroke="$(c)" stroke-width="0.3" transform="matrix($(svgtf))">
     <line x1="0" y1="1" x2="0" y2="-1" />
     <line x1="1" y1="2" x2="1" y2="-2" />
     <line x1="2" y1="3" x2="2" y2="-3" />
     <line x1="3" y1="4" x2="3" y2="-4" />
     <line x1="4" y1="5" x2="4" y2="-5" />
     <line x1="5" y1="6" x2="5" y2="-6" />
     <line x1="6" y1="7" x2="6" y2="-7" />
     <line x1="7" y1="8" x2="7" y2="-8" />
    </g>
    """
    push!(canvas.elements, elem)
    return
end

function drawdashedwedge!(canvas::SvgCanvas, seg, ucolor, vcolor)
    """ u ◁ v """
    ucolor == vcolor && return drawdashedwedge!(canvas, seg, ucolor)
    d = Geometry.distance(seg)
    scalef = Point2D(d / 7, d / 14 * canvas.wedgewidthf)
    rotatef = unitvector(seg)
    translf = seg.u
    svgtf = svgtransform(transformmatrix(scalef, rotatef, translf))
    uc = svgcolor(ucolor)
    vc = svgcolor(vcolor)
    elem = """<g stroke-width="0.3" transform="matrix($(svgtf))">
     <line x1="0" y1="1" x2="0" y2="-1" stroke="$(uc)" />
     <line x1="1" y1="2" x2="1" y2="-2" stroke="$(uc)" />
     <line x1="2" y1="3" x2="2" y2="-3" stroke="$(uc)" />
     <line x1="3" y1="4" x2="3" y2="-4" stroke="$(uc)" />
     <line x1="4" y1="5" x2="4" y2="-5" stroke="$(vc)" />
     <line x1="5" y1="6" x2="5" y2="-6" stroke="$(vc)" />
     <line x1="6" y1="7" x2="6" y2="-7" stroke="$(vc)" />
     <line x1="7" y1="8" x2="7" y2="-8" stroke="$(vc)" />
    </g>
    """
    push!(canvas.elements, elem)
    return
end


function drawwave!(canvas::SvgCanvas, seg, color)
    d = Geometry.distance(seg)
    scalef = Point2D(d / 7, d / 2 * canvas.wedgewidthf)
    rotatef = unitvector(seg)
    translf = seg.u
    svgtf = svgtransform(transformmatrix(scalef, rotatef, translf))
    c = svgcolor(color)
    elem = """<polyline points="0,0 0.5,0 1,1 2,-1 3,1 4,-1 5,1 6,-1 6.5,0 7,0"
     stroke="$(c)" stroke-width="0.2" fill="none" transform="matrix($(svgtf))"/>
    """
    push!(canvas.elements, elem)
    return
end

function drawwave!(canvas::SvgCanvas, seg, ucolor, vcolor)
    ucolor == vcolor && return drawwave!(canvas, seg, ucolor)
    d = Geometry.distance(seg)
    scalef = Point2D(d / 7, d / 2 * canvas.wedgewidthf)
    rotatef = unitvector(seg)
    translf = seg.u
    svgtf = svgtransform(transformmatrix(scalef, rotatef, translf))
    uc = svgcolor(ucolor)
    vc = svgcolor(vcolor)
    elem = """<g stroke-width="0.2" fill="none" transform="matrix($(svgtf))">
     <polyline points="0,0 0.5,0 1,1 2,-1 3,1 3.5,0" stroke="$(uc)" />
     <polyline points="3.5,0 4,-1 5,1 6,-1 6.5,0 7,0" stroke="$(vc)" />
     </g>
    """
    push!(canvas.elements, elem)
    return
end

"""
write svg string to file.
"""
function savesvg(the_svg_string::String, filename::String)
    if ! contains(filename, ".svg")
        filename *= ".svg"
    end

    f = open(filename, "w")
    write(f, the_svg_string)
    close(f)
    return filename
end
