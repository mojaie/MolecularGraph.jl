#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

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


function drawline!(canvas::Canvas, uv::SMatrix{2,2}, color::Color, dashed::Bool)
    option = dashed ? """ stroke-dasharray="10,10" """ : " "
    uvf = crdf.(uv)
    elem = """<line x1="$(uvf[1])" y1="$(uvf[3])"
     x2="$(uvf[2])" y2="$(uvf[4])" stroke="$(rgbf(color))"$(option)/>
    """
    push!(canvas.elements, elem)
    return
end

function drawline!(canvas::Canvas, uv::SMatrix{2,2}, ucolor::Color,
                   vcolor::Color, dashed::Bool)
    option = dashed ? """ stroke-dasharray="10,10" """ : " "
    res = []
    if ucolor == vcolor
        drawline!(canvas, uv, ucolor, dashed)
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

drawline!(canvas, uv, ucolor, vcolor) = drawline!(
    canvas, uv, ucolor, vcolor, false)
drawdashedline!(canvas, uv, ucolor, vcolor) = drawline!(
    canvas, uv, ucolor, vcolor, true)


function drawwedge!(canvas::Canvas, uv::SMatrix{2,2}, color::Color)
    """ u ▶ v """
    uvlen = norm(utov(uv))
    tf = transformmatrix(
        SVector{2}(uvlen, uvlen / 2 * canvas.wedgewidthf),
        normalize(utov(uv)), vecU(uv)
    )
    tfarr = reshape(tf[1:2, :], (1, 6))
    tftxt = join(crdf.(tfarr), " ")
    elem = """<polygon points="0,1 0,-1 1,0" fill="$(rgbf(color))"
     transform="matrix($(tftxt))"/>
    """
    push!(canvas.elements, elem)
    return
end


function drawdashedwedge!(canvas::Canvas, uv::SMatrix{2,2}, color::Color)
    """ u ▶ v """
    uvlen = norm(utov(uv))
    tf = transformmatrix(
        SVector{2}(uvlen / 7, uvlen / 14 * canvas.wedgewidthf),
        normalize(utov(uv)), vecU(uv)
    )
    tfarr = reshape(tf[1:2, :], (1, 6))
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


function drawwave!(canvas::Canvas, uv::SMatrix{2,2}, color::Color)
    uvlen = norm(utov(uv))
    tf = transformmatrix(
        SVector{2}(uvlen / 7, uvlen / 2 * canvas.wavewidthf),
        normalize(utov(uv)), vecU(uv)
    )
    tfarr = reshape(tf[1:2, :], (1, 6))
    tftxt = join(crdf.(tfarr), " ")
    elem = """<polyline points="0,0 0.5,0 1,1 2,-1 3,1 4,-1 5,1 6,-1 6.5,0 7,0"
     stroke="$(rgbf(color))" stroke-width="0.2"
     fill="none" transform="matrix($(tftxt))"/>
    """
    push!(canvas.elements, elem)
    return
end


function drawtext!(canvas::Canvas, pos::SVector{2}, atom::Atom, align=:center)
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
    alg = pos + [xoffset[align], canvas.fontsize / 2]
    algf = crdf.(alg)
    color = Color(getcolor(atom)...)
    svgtext = markup(
        atom, direction[align],
        """<tspan baseline-shift="-25%" font-size="$(small)">""", "</tspan>",
        """<tspan baseline-shift="50%" font-size="$(small)">""", "</tspan>"
    )
    elem = """<text x="$(algf[1])" y="$(algf[2])" font-size="$(canvas.fontsize)"
     fill="$(rgbf(color))"$(option[align])>$(svgtext)</text>
    """
    push!(canvas.elements, elem)
    return
end
