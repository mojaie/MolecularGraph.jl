#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    SvgCanvas,
    tosvg,
    drawsvg,
    initcanvas!


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

    elements::Vector{String}
    coords::Cartesian2D
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

        canvas.elements = []
        canvas.valid = false
        return canvas
    end
end


function tosvg(canvas::Canvas, width::Int, height::Int)
    vbWf = fmt(".2f", canvas.viewboxW)
    vbHf = fmt(".2f", canvas.viewboxH)
    bgc = format("rgb({:d}, {:d}, {:d})", canvas.bgcolor)
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
     fill="$(bgc)" opacity="$(canvas.opacity)"/>
    """
    endsvg = "</svg>\n"
    join([header, bg, canvas.elements..., endsvg], "")
end


"""
    drawsvg(mol::VectorMol, width::Int, height::Int)

Generate molecular structure image as a SVG format string.

`width` and `height` specifies the size of the image (width and height
attribute of svg tag).
"""
function drawsvg(mol::VectorMol, width::Int, height::Int)
    canvas = SvgCanvas()
    draw2d!(canvas, mol)
    tosvg(canvas, width, height)
end


"""
    boundary(mol::VectorMol) -> (top, left, width, height, unit)

Get boundaries and an appropriate bond length unit for the molecule drawing
canvas.
"""
function boundary(mol::VectorMol)
    coords = mol[:coords2d]
    (left, right) = extrema(x_components(coords))
    (bottom, top) = extrema(y_components(coords))
    width = right - left
    height = top - bottom
    dists = []
    # Size unit
    for bond in edgevalues(mol)
        u = _point(coords, bond.u)
        v =  _point(coords, bond.v)
        d = norm(v - u)
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


"""
    initcanvas!(canvas::Canvas, mol::VectorMol)

Move and adjust the size of the molecule for drawing.
"""
function initcanvas!(canvas::Canvas, mol::VectorMol)
    nodecount(mol) == 0 && return
    coords = rawdata(mol[:coords2d])
    (top, left, width, height, unit) = boundary(mol)
    sf = canvas.scalef / unit
    canvas.coords = cartesian2d(
        (coords .- [left top]) .* [1 -1] .* sf
        .+ [canvas.paddingX canvas.paddingY]
    )
    canvas.viewboxW = width * sf + canvas.paddingX * 2
    canvas.viewboxH = height * sf + canvas.paddingY * 2
    canvas.valid = true
    return
end


function trimbond(canvas, seg, uvis, vvis)
    if uvis && vvis
        return trim_uv(seg, canvas.trimoverlapf * 2)
    elseif uvis
        return trim_u(seg, canvas.trimoverlapf)
    elseif vvis
        return trim_v(seg, canvas.trimoverlapf)
    end
    return seg
end


function singlebond!(canvas, coords, ucolor, vcolor, uvis, vvis)
    seg = trimbond(canvas, coords, uvis, vvis)
    drawline!(canvas, seg, ucolor, vcolor)
    return
end


function wedged!(canvas, coords, ucolor, vcolor, uvis, vvis)
    seg = trimbond(canvas, coords, uvis, vvis)
    drawwedge!(canvas, seg, ucolor, vcolor)
    return
end


function dashedwedged!(canvas, coords, ucolor, vcolor, uvis, vvis)
    seg = trimbond(canvas, coords, uvis, vvis)
    drawdashedwedge!(canvas, seg, ucolor, vcolor)
    return
end


function wavesingle!(canvas, coords, ucolor, vcolor, uvis, vvis)
    seg = trimbond(canvas, coords, uvis, vvis)
    drawwave!(canvas, seg, ucolor, vcolor)
    return
end


function doublebond!(canvas, coords, ucolor, vcolor, uvis, vvis)
    dist = canvas.scalef * canvas.mbwidthf / 2
    seg = trimbond(canvas, coords, uvis, vvis)
    seg1 = translate(seg, pi / 2, dist)
    seg2 = translate(seg, -pi / 2, dist)
    drawline!(canvas, seg1, ucolor, vcolor)
    drawline!(canvas, seg2, ucolor, vcolor)
    return
end


function crossdouble!(canvas, coords, ucolor, vcolor, uvis, vvis)
    dist = canvas.scalef * canvas.mbwidthf / 2
    seg = trimbond(canvas, coords, uvis, vvis)
    seg1 = translate(seg, pi / 2, dist)
    seg2 = translate(seg, -pi / 2, dist)
    buf = point(seg1, 2)
    setcoord!(seg1, point(seg2, 2), 2)
    setcoord!(seg2, buf, 2)
    drawline!(canvas, seg1, ucolor, vcolor)
    drawline!(canvas, seg2, ucolor, vcolor)
    return
end


function ringdouble!(canvas, coords, ucolor, vcolor, uvis, vvis, direction)
    seg = trimbond(canvas, coords, uvis, vvis)
    dist = canvas.scalef * canvas.mbwidthf
    segin = trim_uv(translate(seg, direction, dist), canvas.triminnerf)
    drawline!(canvas, seg, ucolor, vcolor)
    drawline!(canvas, segin, ucolor, vcolor)
    return
end

# NOTE: the direction is reversed due to x-axis reflection
clockwisedouble!(canvas, coords, ucolor, vcolor, uvis, vvis
    ) = ringdouble!(canvas, coords, ucolor, vcolor, uvis, vvis, pi / 2)
counterdouble!(canvas, coords, ucolor, vcolor, uvis, vvis
    ) = ringdouble!(canvas, coords, ucolor, vcolor, uvis, vvis, -pi / 2)


function triplebond!(canvas, coords, ucolor, vcolor, uvis, vvis)
    dist = canvas.scalef * canvas.mbwidthf
    seg = trimbond(canvas, coords, uvis, vvis)
    seg1 = translate(seg, pi / 2, dist)
    seg2 = translate(seg, -pi / 2, dist)
    drawline!(canvas, seg, ucolor, vcolor)
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
    canvas, order, notation, coords, ucolor, vcolor, uvis, vvis
) = BOND_DRAWER[order][notation](canvas, coords, ucolor, vcolor, uvis, vvis)


function atomsymbol!(canvas, coords, atomsymbol, color, implicith, charge;
                     direction=:right, anchor=""" text-anchor="end" """,
                     xoffset=0)
    pos = coords + [xoffset, canvas.fontsize/2]
    (x, y) = fmt.(".2f", pos)
    small = round(Int, canvas.fontsize * 0.7)
    content = atommarkup(
        atomsymbol, charge, implicith, direction,
        """<tspan baseline-shift="-25%" font-size="$(small)">""", "</tspan>",
        """<tspan baseline-shift="50%" font-size="$(small)">""", "</tspan>"
    )
    c = format("rgb({:d}, {:d}, {:d})", color)
    elem = """<text x="$(x)" y="$(y)" font-size="$(canvas.fontsize)"
     fill="$(c)"$(anchor)>$(content)</text>
    """
    push!(canvas.elements, elem)
    return
end

setatomright!(
    canvas, coords, sym, color, hcnt, charge
) = atomsymbol!(
    canvas, coords, sym, color, hcnt, charge,
    direction=:left, xoffset=canvas.fontsize/2
)

setatomcenter!(
    canvas, coords, sym, color, hcnt, charge
) = atomsymbol!(
    canvas, coords, sym, color, hcnt, charge,
    anchor=""" text-anchor="middle" """
)

setatomleft!(
    canvas, coords, sym, color, hcnt, charge
) = atomsymbol!(
    canvas, coords, sym, color, hcnt, charge,
    anchor=" ", xoffset=canvas.fontsize/-2
)


function setatomnote!(canvas, coords, text, color, bgcolor)
    size = round(Int, canvas.fontsize * canvas.annotsizef)
    (bx, by) = fmt.(".2f", coords)
    (tx, ty) = fmt.(".2f", coords + [0 size])
    c = format("rgb({:d}, {:d}, {:d})", color)
    bc = format("rgb({:d}, {:d}, {:d})", bgcolor)
    elem = """<g>
     <rect x="$(bx)" y="$(by)" width="$(size)" height="$(size)" rx="$(size/2)" ry="$(size/2)" fill="$(bc)" />
     <text x="$(tx)" y="$(ty)" font-size="$(size)" fill="$(c)">$(text)</text>
    </g>
    """
    push!(canvas.elements, elem)
    return
end



function drawline!(canvas, seg, color; isdashed=false)
    option = isdashed ? """ stroke-dasharray="10,10" """ : " "
    ((ux, uy), (vx, vy)) = fmt(".2f", seg)
    c = format("rgb({:d}, {:d}, {:d})", color)
    elem = """<line x1="$(ux)" y1="$(uy)" x2="$(vx)" y2="$(vy)" stroke="$(c)"$(option)/>"""
    push!(canvas.elements, elem)
    return
end

function drawline!(canvas, seg, ucolor, vcolor; isdashed=false)
    if ucolor == vcolor
        drawline!(canvas, seg, ucolor, isdashed=isdashed)
        return
    end
    option = isdashed ? """ stroke-dasharray="10,10" """ : " "
    ((ux, uy), (vx, vy)) = fmt(".2f", seg)
    (mx, my) = fmt(".2f", midpoint(seg))
    uc = format("rgb({:d}, {:d}, {:d})", ucolor)
    vc = format("rgb({:d}, {:d}, {:d})", vcolor)
    elem = """<line x1="$(ux)" y1="$(uy)" x2="$(mx)" y2="$(my)" stroke="$(uc)"$(option)/>
    <line x1="$(mx)" y1="$(my)" x2="$(vx)" y2="$(vy)" stroke="$(vc)"$(option)/>
    """
    push!(canvas.elements, elem)
    return
end

drawdashedline!(canvas, seg, ucolor, vcolor) = drawline!(
    canvas, seg, ucolor, vcolor, isdashed=true)


function drawwedge!(canvas, seg, color)
    """ u ◀︎ v """
    uvlen = norm(_vector(seg))
    tf = transformmatrix(
        point(uvlen, uvlen / 2 * canvas.wedgewidthf),
        point(normalize(_vector(seg))...),
        point(_u(seg)...)
    )
    mat = join(fmt.(".2f", tf[[1, 2, 4, 5, 7, 8]]), " ")
    c = format("rgb({:d}, {:d}, {:d})", color)
    elem = """<polygon points="0,0 1,1 1,-1" fill="$(c)" transform="matrix($(mat))"/>
    """
    push!(canvas.elements, elem)
    return
end

function drawwedge!(canvas, seg, ucolor, vcolor)
    """ u ◀︎ v """
    if ucolor == vcolor
        drawwedge!(canvas, seg, ucolor)
        return
    end
    uvlen = norm(_vector(seg))
    tf = transformmatrix(
        point(uvlen, uvlen / 2 * canvas.wedgewidthf),
        point(normalize(_vector(seg))...),
        point(_u(seg)...)
    )
    mat = join(fmt.(".2f", tf[[1, 2, 4, 5, 7, 8]]), " ")
    uc = format("rgb({:d}, {:d}, {:d})", ucolor)
    vc = format("rgb({:d}, {:d}, {:d})", vcolor)
    elem = """<g stroke-width="0.3" transform="matrix($(mat))">
     <polygon points="0,0 0.5,-0.5 0.5,0.5" fill="$(uc)"/>
     <polygon points="0.5,-0.5 0.5,0.5 1,1 1,-1" fill="$(vc)"/>
     </g>
    """
    push!(canvas.elements, elem)
    return
end


function drawdashedwedge!(canvas, seg, color)
    """ u ◁ v """
    uvlen = norm(_vector(seg))
    tf = transformmatrix(
        point(uvlen / 7, uvlen / 14 * canvas.wedgewidthf),
        point(normalize(_vector(seg))...),
        point(_u(seg)...)
    )
    mat = join(fmt.(".2f", tf[[1, 2, 4, 5, 7, 8]]), " ")
    c = format("rgb({:d}, {:d}, {:d})", color)
    elem = """<g stroke="$(c)" stroke-width="0.3" transform="matrix($(mat))">
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

function drawdashedwedge!(canvas, seg, ucolor, vcolor)
    """ u ◁ v """
    if ucolor == vcolor
        drawdashedwedge!(canvas, seg, ucolor)
        return
    end
    uvlen = norm(_vector(seg))
    tf = transformmatrix(
        point(uvlen / 7, uvlen / 14 * canvas.wedgewidthf),
        point(normalize(_vector(seg))...),
        point(_u(seg)...)
    )
    mat = join(fmt.(".2f", tf[[1, 2, 4, 5, 7, 8]]), " ")
    uc = format("rgb({:d}, {:d}, {:d})", ucolor)
    vc = format("rgb({:d}, {:d}, {:d})", vcolor)
    elem = """<g stroke-width="0.3" transform="matrix($(mat))">
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


function drawwave!(canvas, seg, color)
    uvlen = norm(_vector(seg))
    tf = transformmatrix(
        point(uvlen / 7, uvlen / 2 * canvas.wedgewidthf),
        point(normalize(_vector(seg))...),
        point(_u(seg)...)
    )
    mat = join(fmt.(".2f", tf[[1, 2, 4, 5, 7, 8]]), " ")
    c = format("rgb({:d}, {:d}, {:d})", color)
    elem = """<polyline points="0,0 0.5,0 1,1 2,-1 3,1 4,-1 5,1 6,-1 6.5,0 7,0"
     stroke="$(c)" stroke-width="0.2" fill="none" transform="matrix($(mat))"/>
    """
    push!(canvas.elements, elem)
    return
end

function drawwave!(canvas, seg, ucolor, vcolor)
    if ucolor == vcolor
        drawwave!(canvas, seg, ucolor)
        return
    end
    uvlen = norm(_vector(seg))
    tf = transformmatrix(
        point(uvlen / 7, uvlen / 2 * canvas.wedgewidthf),
        point(normalize(_vector(seg))...),
        point(_u(seg)...)
    )
    mat = join(fmt.(".2f", tf[[1, 2, 4, 5, 7, 8]]), " ")
    uc = format("rgb({:d}, {:d}, {:d})", ucolor)
    vc = format("rgb({:d}, {:d}, {:d})", vcolor)
    elem = """<g stroke-width="0.2" fill="none" transform="matrix($(mat))">
     <polyline points="0,0 0.5,0 1,1 2,-1 3,1 3.5,0" stroke="$(uc)" />
     <polyline points="3.5,0 4,-1 5,1 6,-1 6.5,0 7,0" stroke="$(vc)" />
     </g>
    """
    push!(canvas.elements, elem)
    return
end
