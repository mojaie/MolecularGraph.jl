#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    SvgCanvas,
    tosvg,
    drawsvg!,
    initcanvas!,
    drawline!,
    drawwedge!,
    drawdashedwedge!,
    drawwave!,
    drawsymbol!,
    drawsymbolright!,
    drawsymbolcenter!,
    drawsymbolleft!,
    drawtext!


rgbf(c::Color) = @sprintf "rgb(%d, %d, %d)" c.r c.g c.b
crdf(c::Real) = @sprintf "%.2f" c


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
    scalef::Float64
    paddingX::Float64
    paddingY::Float64

    coords::Matrix{Float64}
    viewboxW::Float64
    viewboxH::Float64

    elements::Vector{String}
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
        canvas.scalef = 30.0 # suitable for fontsize=14 and default line width
        canvas.paddingX = 30.0
        canvas.paddingY = 30.0

        canvas.elements = []
        canvas.valid = false
        canvas
    end
end


function tosvg(canvas::Canvas, width, height)
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


function drawsvg!(mol::VectorMol, width, height)
    draw2d_annot!(mol)
    canvas = SvgCanvas()
    draw2d!(canvas, mol)
    tosvg(canvas, width, height)
end


function initcanvas!(canvas::Canvas, mol::VectorMol)
    if atomcount(mol) == 0
        return
    end
    coords = mol[:Coords2D]
    (top, left, width, height, unit) = boundary(mol, coords)
    sf = canvas.scalef / unit
    coords = (coords .- [left top]) .* [1 -1] .* sf
    canvas.coords = coords .+ [canvas.paddingX canvas.paddingY]
    canvas.viewboxW = width * sf + canvas.paddingX * 2
    canvas.viewboxH = height * sf + canvas.paddingY * 2
    canvas.valid = true
    return
end


function drawline!(canvas, uv, color::Color, isdashed::Bool)
    option = isdashed ? """ stroke-dasharray="10,10" """ : " "
    uvf = crdf.(uv)
    elem = """<line x1="$(uvf[1])" y1="$(uvf[3])"
     x2="$(uvf[2])" y2="$(uvf[4])" stroke="$(rgbf(color))"$(option)/>
    """
    push!(canvas.elements, elem)
    return
end

function drawline!(canvas, uv, ucolor, vcolor, isdashed)
    option = isdashed ? """ stroke-dasharray="10,10" """ : " "
    res = []
    if ucolor == vcolor
        drawline!(canvas, uv, ucolor, isdashed)
        return
    end
    uvf = crdf.(uv)
    midf = crdf.((vecU(uv) + vecV(uv)) / 2)
    elem = """<line x1="$(uvf[1])" y1="$(uvf[3])"
     x2="$(midf[1])" y2="$(midf[2])" stroke="$(rgbf(ucolor))"$(option)/>
    <line x1="$(midf[1])" y1="$(midf[2])"
     x2="$(uvf[2])" y2="$(uvf[4])" stroke="$(rgbf(vcolor))"$(option)/>
    """
    push!(canvas.elements, elem)
    return
end

drawline!(canvas, uv, ucolor::Color, vcolor::Color) = drawline!(
    canvas, uv, ucolor, vcolor, false)
drawdashedline!(canvas, uv, ucolor, vcolor) = drawline!(
    canvas, uv, ucolor, vcolor, true)


function drawwedge!(canvas, uv, color)
    """ u ▶ v """
    uvlen = norm(utov(uv))
    tf = transformmatrix(
        SVector{2}(uvlen, uvlen / 2 * canvas.wedgewidthf),
        normalize(utov(uv)), vecU(uv)
    )
    # TODO: StaticArrays bug?
    tfarr = tf[[1, 2, 4, 5, 7, 8]]
    # tfarr = reshape(tf[1:2, :], (1, 6))
    tftxt = join(crdf.(tfarr), " ")
    elem = """<polygon points="0,1 0,-1 1,0" fill="$(rgbf(color))"
     transform="matrix($(tftxt))"/>
    """
    push!(canvas.elements, elem)
    return
end


function drawdashedwedge!(canvas, uv, color)
    """ u ▶ v """
    uvlen = norm(utov(uv))
    tf = transformmatrix(
        SVector{2}(uvlen / 7, uvlen / 14 * canvas.wedgewidthf),
        normalize(utov(uv)), vecU(uv)
    )
    # TODO: StaticArrays bug?
    tfarr = tf[[1, 2, 4, 5, 7, 8]]
    # tfarr = reshape(tf[1:2, :], (1, 6))
    tftxt = join(crdf.(tfarr), " ")
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


function drawwave!(canvas, uv, color)
    uvlen = norm(utov(uv))
    tf = transformmatrix(
        SVector{2}(uvlen / 7, uvlen / 2 * canvas.wavewidthf),
        normalize(utov(uv)), vecU(uv)
    )
    # TODO: StaticArrays bug?
    tfarr = tf[[1, 2, 4, 5, 7, 8]]
    # tfarr = reshape(tf[1:2, :], (1, 6))
    tftxt = join(crdf.(tfarr), " ")
    elem = """<polyline points="0,0 0.5,0 1,1 2,-1 3,1 4,-1 5,1 6,-1 6.5,0 7,0"
     stroke="$(rgbf(color))" stroke-width="0.2"
     fill="none" transform="matrix($(tftxt))"/>
    """
    push!(canvas.elements, elem)
    return
end


function drawsymbol!(canvas, pos, symbol, charge, hcount, color,
                     direction, option, xoffset)
    tpos = crdf.(pos + [xoffset, canvas.fontsize / 2])
    small = round(Int, canvas.fontsize * 0.7)
    atom = atomsvg(symbol, charge, hcount, direction, small)
    elem = """<text x="$(tpos[1])" y="$(tpos[2])" font-size="$(canvas.fontsize)"
     fill="$(rgbf(color))"$(option)>$(atom)</text>
    """
    push!(canvas.elements, elem)
    return
end

function drawsymbolright!(canvas, pos, symbol, charge, hcount, color)
    drawsymbol!(canvas, pos, symbol, charge, hcount, color,
                :left, """ text-anchor="end" """, canvas.fontsize / 2)
end

function drawsymbolcenter!(canvas, pos, symbol, charge, hcount, color)
    drawsymbol!(canvas, pos, symbol, charge, hcount, color,
                :right, """ text-anchor="middle" """, 0)
end

function drawsymbolleft!(canvas, pos, symbol, charge, hcount, color)
    drawsymbol!(canvas, pos, symbol, charge, hcount, color,
                :right, " ", canvas.fontsize / -2)
end


function drawtext!(canvas, pos, text, fontsizef, color, bgcolor)
    size = round(Int, canvas.fontsize * fontsizef)
    rpos = crdf.(pos)
    tpos = crdf.(pos + [0, size])
    elem = """<g>
        <rect x="$(rpos[1])" y="$(rpos[2])" width="$(size)" height="$(size)"
         rx="$(size/2)" ry="$(size/2)" fill="$(rgbf(bgcolor))" />
        <text x="$(tpos[1])" y="$(tpos[2])" font-size="$(size)"
         fill="$(rgbf(color))">$(text)</text>
    </g>
    """
    push!(canvas.elements, elem)
    return
end

function atomsvg(symbol, charge, hcount, direction, small)
    atommarkup(
        symbol, charge, hcount, direction,
        """<tspan baseline-shift="-25%" font-size="$(small)">""", "</tspan>",
        """<tspan baseline-shift="50%" font-size="$(small)">""", "</tspan>"
    )
end
