#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

using Printf

export
    SvgCanvas,
    tosvg,
    drawsvg!,
    initialize!,
    drawline!,
    drawwedge!,
    drawdashedwedge!,
    drawwave!,
    drawtext!


mutable struct SvgCanvas <: Canvas
    fontweight::String
    fontfamily::String
    fontsize::UInt8
    bgcolor::Color
    opacity::UInt8

    mbwidthf::Float32
    wedgewidthf::Float32
    wavewidthf::Float32
    triminnerf::Float32
    trimoverlapf::Float32

    coords::Matrix{Float32}
    viewboxW::Float32
    viewboxH::Float32

    elements::Vector{String}

    function SvgCanvas()
        canvas = new()
        canvas.fontweight = "normal"
        canvas.fontfamily = "Helvetica"
        canvas.fontsize = 14
        canvas.bgcolor = Color(255, 255, 255)
        canvas.opacity = 0

        canvas.mbwidthf = 0.15
        canvas.wedgewidthf = 0.2
        canvas.wavewidthf = 0.2
        canvas.triminnerf = 0.2
        canvas.trimoverlapf = 0.3

        canvas.elements = []
        canvas
    end
end


function tosvg(canvas::Canvas, width::Real, height::Real)
    header = """<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:xlink="http://www.w3.org/1999/xlink"
     version="1.2" baseProfile="tiny"
     text-rendering="geometricPrecision"
     preserveAspectRatio="xMidYMid meet"
     font-weight="$(canvas.fontweight)"
     font-family="$(canvas.fontfamily)"
     width="$(width)" height="$(height)"
     viewBox="0 0 $(canvas.viewboxW) $(canvas.viewboxH)">\n
    """
    bg = """<rect x="0", y="0"
     width="$(canvas.viewboxW)" height="$(canvas.viewboxH)"
     fill="$(canvas.bgcolor)" opacity="$(canvas.opacity)"/>\n
    """
    endsvg = "</svg>"
    join([header, bg, canvas.elements..., endsvg], "")
end


function drawsvg!(mol::MolecularGraph, width::Real, height::Real)
    readytodraw!(mol)
    canvas = SvgCanvas()
    draw!(canvas, mol)
    tosvg(canvas, width, height)
end


function initialize!(canvas::Canvas, mol::MolecularGraph)
    coords = coordsvector(mol)
    (top, left, width, height) = boundary(coords)
    # suitable for fontsize=14 and default line width
    sizefactor = 30 / sizeunit(mol, coords)
    bottom = top - height
    canvas.coords = coords .- [left bottom] .* [1 -1] .* sizefactor
    canvas.viewboxW = width * sizefactor
    canvas.viewboxH = height * sizefactor
    return
end


function drawline!(canvas::Canvas, seg::Segment, color::Color, dashed::Bool)
    option = dashed ? """ stroke-dasharray="10,10" """ : " "
    u, v = (seg.u, seg.v)
    elem = """<line x1="$(u.x)" y1="$(u.y)" x2="$(v.x)" y2="$(v.y)"
     stroke="$(color)"$(option)/>\n
    """
    push!(canvas.elements, elem)
    return
end

function drawline!(canvas::Canvas, seg::Segment, ucolor::Color,
                   vcolor::Color, dashed::Bool)
    option = dashed ? """ stroke-dasharray="10,10" """ : " "
    u, v = (seg.u, seg.v)
    res = []
    if ucolor == vcolor
        drawline!(canvas, seg, ucolor, dashed)
        return
    end
    mid = (u + v) / 2
    elem = """<line x1="$(u.x)" y1="$(u.y)" x2="$(mid.x)" y2="$(mid.y)"
     stroke="$(ucolor)"$(option)/>\n
    <line x1="$(mid.x)" y1="$(mid.y)" x2="$(v.x)" y2="$(v.y)"
     stroke="$(vcolor)"$(option)/>\n
    """
    push!(canvas.elements, elem)
    return
end

drawline!(canvas, seg, ucolor, vcolor) = drawline!(
    canvas, seg, ucolor, vcolor, false)
drawdashedline!(canvas, seg, ucolor, vcolor) = drawline!(
    canvas, seg, ucolor, vcolor, true)


function drawwedge!(canvas::Canvas, seg::Segment, color::Color)
    """ u ▶ v """
    seglen = length(seg)
    tf = transformmatrix(
        seglen, seglen / 2 * canvas.mbwidthf,
        (seg.v.x - seg.u.x) / seglen, (seg.v.y - seg.u.y) / seglen,
        seg.u.x, seg.v.y
    )
    tfarr = reshape(tf, (1, 6))
    tftxt = join([@sprintf("%.2f", v) for v in tfarr], " ")
    elem = """<polygon points="0,1 0,-1, 1,0" fill="$(color)"
     transform="matrix($(tftxt))"/>\n
    """
    push!(canvas.elements, elem)
    return
end


function drawdashedwedge!(canvas::Canvas, seg::Segment, color::Color)
    """ u ▶ v """
    seglen = length(seg)
    tf = transformmatrix(
        seglen / 7, seglen / 7 * canvas.wavewidthf,
        (seg.v.x - seg.u.x) / seglen, (seg.v.y - seg.u.y) / seglen,
        seg.u.x, seg.v.y
    )
    tfarr = reshape(tf, (1, 6))
    tftxt = join([@sprintf("%.2f", v) for v in tfarr], " ")
    elem = """<g stroke="$(color)" transform="matrix($(tftxt))">\n
     <line x1="0" y1="7" x2="0" y2="-7" />\n
     <line x1="1" y1="6" x2="0" y2="-6" />\n
     <line x1="2" y1="5" x2="0" y2="-5" />\n
     <line x1="3" y1="4" x2="0" y2="-4" />\n
     <line x1="4" y1="3" x2="0" y2="-3" />\n
     <line x1="5" y1="2" x2="0" y2="-2" />\n
     <line x1="6" y1="1" x2="0" y2="-1" />\n
    </g>
    """
    push!(canvas.elements, elem)
    return
end


function drawwave!(canvas::Canvas, seg::Segment, color::Color)
    seglen = length(seg)
    tf = transformmatrix(
        seglen / 7, seglen / 2 * canvas.wavewidthf,
        (seg.v.x - seg.u.x) / seglen, (seg.v.y - seg.u.y) / seglen,
        seg.u.x, seg.v.y
    )
    tfarr = reshape(tf, (1, 6))
    tftxt = join([@sprintf("%.2f", v) for v in tfarr], " ")
    elem = """<polyline points="0,0 1,1, 2,-1 3,1 4,-1 5,1 6,-1 7,0"
     stroke="$(color)" fill="none" transform="matrix($(tftxt))"/>\n
    """
    push!(canvas.elements, elem)
    return
end


function drawtext!(canvas::Canvas, pos::Point2D, text::String, color::Color,
                   align=:center)
    small = canvas.fontsize * 0.7
    svgtxt = text[:]
    replace(svgtxt,
        r"<sub>" => """<tspan baseline-shift="-25%" font-size="$(small)">""")
    replace(svgtxt,
        r"<sup>" => """<tspan baseline-shift="50%" font-size="$(small)">""")
    replace(svgtxt, r"</sub>" => "</tspan>")
    replace(svgtxt, r"</sup>" => "</tspan>")
    option = Dict(
        :center => """ text-anchor="middle" """,
        :right => """ text-anchor="end" """, :left => " ",
    )
    xoffset = Dict(
        :center => 0, :right => canvas.fontsize / 2,
        :left => -canvas.fontsize / 2
    )
    px = pos.x + xoffset[align]
    py = pos.y + canvas.fontsize / 2
    elem = """<text x="$(px)" y="$(py)" font-size="$(canvas.fontsize)"
     fill="$(color)"$(option)>$(svgtxt)</text>\n
    """
    push!(canvas.elements, elem)
    return
end
