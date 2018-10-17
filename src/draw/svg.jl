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


rgbf(c::Color) = @sprintf "rgb(%d, %d, %d)" c.r c.g c.b
crdf(c::Real) = @sprintf "%.2f" c


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
    scalef::Float32
    paddingX::Float32
    paddingY::Float32

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
        canvas.wedgewidthf = 0.3
        canvas.wavewidthf = 0.2
        canvas.triminnerf = 0.2
        canvas.trimoverlapf = 0.3
        canvas.scalef = 30 # suitable for fontsize=14 and default line width
        canvas.paddingX = 30
        canvas.paddingY = 30

        canvas.elements = []
        canvas
    end
end


function tosvg(canvas::Canvas, width::Real, height::Real)
    vbWf = crdf(canvas.viewboxW)
    vbHf = crdf(canvas.viewboxH)
    header = """<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:xlink="http://www.w3.org/1999/xlink"
     version="1.2" baseProfile="tiny"
     text-rendering="geometricPrecision"
     preserveAspectRatio="xMidYMid meet"
     font-weight="$(canvas.fontweight)"
     font-family="$(canvas.fontfamily)"
     width="$(width)" height="$(height)"
     viewBox="0 0 $(vbWf) $(vbHf)">
    """
    bg = """<rect x="0" y="0" width="$(vbWf)" height="$(vbHf)"
     fill="$(rgbf(canvas.bgcolor))" opacity="$(canvas.opacity)"/>
    """
    endsvg = "</svg>\n"
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
    sf = canvas.scalef / sizeunit(mol, coords)
    coords = (coords .- [left top]) .* [1 -1] .* sf
    canvas.coords = coords .+ [canvas.paddingX canvas.paddingY]
    canvas.viewboxW = width * sf + canvas.paddingX * 2
    canvas.viewboxH = height * sf + canvas.paddingY * 2
    return
end


function drawline!(canvas::Canvas, seg::Segment, color::Color, dashed::Bool)
    option = dashed ? """ stroke-dasharray="10,10" """ : " "
    u, v = (seg.u, seg.v)
    elem = """<line x1="$(crdf(u.x))" y1="$(crdf(u.y))"
     x2="$(crdf(v.x))" y2="$(crdf(v.y))" stroke="$(rgbf(color))"$(option)/>
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
    elem = """<line x1="$(crdf(u.x))" y1="$(crdf(u.y))"
     x2="$(crdf(mid.x))" y2="$(crdf(mid.y))" stroke="$(rgbf(ucolor))"$(option)/>
    <line x1="$(crdf(mid.x))" y1="$(crdf(mid.y))"
     x2="$(crdf(v.x))" y2="$(crdf(v.y))" stroke="$(rgbf(vcolor))"$(option)/>
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
        seglen, seglen / 2 * canvas.wedgewidthf,
        (seg.v.x - seg.u.x) / seglen, (seg.u.y - seg.v.y) / seglen,
        seg.u.x, seg.u.y
    )
    tfarr = reshape(tf, (1, 6))
    tftxt = join([crdf(v) for v in tfarr], " ")
    elem = """<polygon points="0,1 0,-1 1,0" fill="$(rgbf(color))"
     transform="matrix($(tftxt))"/>
    """
    push!(canvas.elements, elem)
    return
end


function drawdashedwedge!(canvas::Canvas, seg::Segment, color::Color)
    """ u ▶ v """
    seglen = length(seg)
    tf = transformmatrix(
        seglen / 7, seglen / 14 * canvas.wedgewidthf,
        (seg.v.x - seg.u.x) / seglen, (seg.u.y - seg.v.y) / seglen,
        seg.u.x, seg.u.y
    )
    tfarr = reshape(tf, (1, 6))
    tftxt = join([crdf(v) for v in tfarr], " ")
    elem = """<g stroke="$(rgbf(color))" stroke-width="0.3" transform="matrix($(tftxt))">
     <line x1="0" y1="8" x2="0" y2="-8" />
     <line x1="1" y1="7" x2="1" y2="-7" />
     <line x1="2" y1="6" x2="2" y2="-6" />
     <line x1="3" y1="5" x2="3" y2="-5" />
     <line x1="4" y1="4" x2="4" y2="-4" />
     <line x1="5" y1="3" x2="5" y2="-3" />
     <line x1="6" y1="2" x2="6" y2="-2" />
     <line x1="7" y1="1" x2="7" y2="-1" />
    </g>
    """
    push!(canvas.elements, elem)
    return
end


function drawwave!(canvas::Canvas, seg::Segment, color::Color)
    seglen = length(seg)
    tf = transformmatrix(
        seglen / 7, seglen / 2 * canvas.wavewidthf,
        (seg.v.x - seg.u.x) / seglen, (seg.u.y - seg.v.y) / seglen,
        seg.u.x, seg.u.y
    )
    tfarr = reshape(tf, (1, 6))
    tftxt = join([crdf(v) for v in tfarr], " ")
    elem = """<polyline points="0,0 0.5,0 1,1 2,-1 3,1 4,-1 5,1 6,-1 6.5,0 7,0"
     stroke="$(rgbf(color))" stroke-width="0.2"
     fill="none" transform="matrix($(tftxt))"/>
    """
    push!(canvas.elements, elem)
    return
end


function drawtext!(canvas::Canvas, pos::Point2D, atom::Atom, align=:center)
    small = round(Int, canvas.fontsize * 0.7)
    direction = Dict(:center => :right, :right => :left, :left => :right)
    option = Dict(
        :center => """ text-anchor="middle" """,
        :right => """ text-anchor="end" """, :left => " ",
    )
    xoffset = Dict(
        :center => 0, :right => canvas.fontsize / 2,
        :left => canvas.fontsize / -2
    )
    px = pos.x + xoffset[align]
    py = pos.y + canvas.fontsize / 2
    color = Color(getcolor(atom)...)
    svgtext = markup(
        atom, direction[align],
        """<tspan baseline-shift="-25%" font-size="$(small)">""", "</tspan>",
        """<tspan baseline-shift="50%" font-size="$(small)">""", "</tspan>"
    )
    elem = """<text x="$(crdf(px))" y="$(crdf(py))" font-size="$(canvas.fontsize)"
     fill="$(rgbf(color))"$(option[align])>$(svgtext)</text>
    """
    push!(canvas.elements, elem)
    return
end
