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


function drawsvg!(mol::MolecularGraph)
    readytodraw!(mol)
    canvas = SvgCanvas()
    draw!(canvas, mol)
    tosvg(canvas)
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


function drawline!(canvas::Canvas, seg::Segment, color::Color;
                  vcolor=nothing, dashed=false)
    option = dashed ? """ stroke-dasharray="10,10" """ : " "
    u, v = (seg.u, seg.v)
    res = []
    if vcolor !== nothing && color != vcolor
        mid = (u + v) / 2
        push!(canvas.elements,
        """<line x1="$(u.x)" y1="$(u.y)" x2="$(mid.x)" y2="$(mid.y)"
         stroke="$(color)"$(option)/>\n
        <line x1="$(mid.x)" y1="$(mid.y)" x2="$(v.x)" y2="$(v.y)"
         stroke="$(vcolor)"$(option)/>\n""")
    else
        push!(canvas.elements,
        """<line x1="$(u.x)" y1="$(u.y)" x2="$(v.x)" y2="$(v.y)"
         stroke="$(color)"$(option)/>\n""")
    end
    return
end

drawline!(canvas, seg, color, vcolor) = drawline!(
    canvas, seg, color, vcolor=vcolor, dashed=false)

drawdashedline!(canvas, seg, color, vcolor) = drawline!(
    canvas, seg, color, vcolor=vcolor, dashed=true)


function drawwedge!(canvas::Canvas, seg::Segment, color::Color)
    """ u ▶ v """
    seglen = length(seg)
    tf = transformmatrix(
        seglen, seglen / 2 * scale.mbwidthf,
        (seg.v.x - seg.u.x) / seglen, (seg.v.y - seg.u.y) / seglen,
        seg.u.x, seg.v.y
    )
    tfarr = reshape(tf, (1, 6))
    tftxt = join([@sprintf("%.2f", v) for v in tfarr], " ")
    push!(canvas.elements,
    """<polygon points="0,1 0,-1, 1,0" fill="$(color)"
     transform="matrix($(tftxt))"/>\n""")
    return
end


function drawdashedwedge!(canvas::Canvas, seg::Segment, color::Color)
    """ u ▶ v """
    seglen = length(seg)
    tf = transformmatrix(
        seglen / 7, seglen / 7 * scale.wavewidthf,
        (seg.v.x - seg.u.x) / seglen, (seg.v.y - seg.u.y) / seglen,
        seg.u.x, seg.v.y
    )
    tfarr = reshape(tf, (1, 6))
    tftxt = join([@sprintf("%.2f", v) for v in tfarr], " ")
    push!(canvas.elements,
    """<g stroke="$(color)" transform="matrix($(tftxt))">\n
     <line x1="0" y1="7" x2="0" y2="-7" />\n
     <line x1="1" y1="6" x2="0" y2="-6" />\n
     <line x1="2" y1="5" x2="0" y2="-5" />\n
     <line x1="3" y1="4" x2="0" y2="-4" />\n
     <line x1="4" y1="3" x2="0" y2="-3" />\n
     <line x1="5" y1="2" x2="0" y2="-2" />\n
     <line x1="6" y1="1" x2="0" y2="-1" />\n
    </g>
    """)
    return
end


function drawwave!(canvas::Canvas, seg::Segment, color::Color)
    seglen = length(seg)
    tf = transformmatrix(
        seglen / 7, seglen / 2 * scale.wavewidthf,
        (seg.v.x - seg.u.x) / seglen, (seg.v.y - seg.u.y) / seglen,
        seg.u.x, seg.v.y
    )
    tfarr = reshape(tf, (1, 6))
    tftxt = join([@sprintf("%.2f", v) for v in tfarr], " ")
    push!(canvas.elements,
    """<polyline points="0,0 1,1, 2,-1 3,1 4,-1 5,1 6,-1 7,0"
     stroke="$(color)" fill="none" transform="matrix($(tftxt))"/>\n""")
    return
end


function drawtext!(canvas::Canvas, pos::Point2D, text::String, color::Color,
                   align=:center)
    small = canvas.fontsize * 0.7
    svgtxt = copy(text)
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
    push!(canvas.elements,
    """<text x="$(px)" y="$(py)" font-size="$(canvas.fontsize)"
     fill="$(color)"$(option)>$(svgtxt)</text>\n""")
    return
end
