#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#


mutable struct SvgCanvas <: Canvas
    scaleunit::Float64
    mbwidthf::Float64
    wedgewidthf::Float64
    wavewidthf::Float64
    triminnerf::Float64
    trimoverlapf::Float64
    linehlwidthf::Float64
    annotsizef::Float64
    hlsizef::Float64
    paddingXf::Float64
    paddingYf::Float64

    fontweight::String
    fontfamily::String
    fontsize::Float64
    bgcolor::RGB
    bgopacity::Float64

    viewboxW::Float64
    viewboxH::Float64
    bgelements::Vector{String}
    elements::Vector{String}
    coords::Matrix{Float64}

    function SvgCanvas(bgcolor::RGB, bgopacity::Float64)
        canvas = new()

        # Geometry
        canvas.scaleunit = 30.0
        canvas.mbwidthf = 0.15
        canvas.wedgewidthf = 0.3
        canvas.wavewidthf = 0.2
        canvas.triminnerf = 0.2
        canvas.trimoverlapf = 0.3
        canvas.linehlwidthf = 0.3
        canvas.annotsizef = 0.7
        canvas.hlsizef = 1.2
        canvas.paddingXf = 1.0
        canvas.paddingYf = 1.0

        # Appearance
        canvas.fontweight = "normal"
        canvas.fontfamily = "Helvetica"
        canvas.fontsize = 14.0
        canvas.bgcolor = bgcolor
        canvas.bgopacity = bgopacity

        # Canvas state
        canvas.viewboxW = 1.0
        canvas.viewboxH = 1.0
        canvas.bgelements = []
        canvas.elements = []

        return canvas
    end
end


svgcoords(p::Point2D) = @sprintf "x=\"%.2f\" y=\"%.2f\"" p.x p.y
svgcirclecoords(p::Point2D) = @sprintf "cx=\"%.2f\" cy=\"%.2f\"" p.x p.y
svgcoords(s::Segment
    ) = @sprintf "x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\"" s.u.x s.u.y s.v.x s.v.y
svgtransform(tf::Array{Float64,2}
    ) = @sprintf "%.2f %.2f %.2f %.2f %.2f %.2f" tf[1] tf[2] tf[4] tf[5] tf[7] tf[8]


function tosvg(canvas::SvgCanvas; width="100%", height="100%", viewbox=true, kwargs...)
    vbWf = @sprintf "%.2f" canvas.viewboxW
    vbHf = @sprintf "%.2f" canvas.viewboxH
    vb = viewbox ? " viewBox=\"0 0 $(vbWf) $(vbHf)\"" : ""
    w = viewbox ? width : vbWf
    h = viewbox ? height : vbHf
    header = """<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:xlink="http://www.w3.org/1999/xlink"
     version="1.1" baseProfile="tiny"
     text-rendering="geometricPrecision"
     preserveAspectRatio="xMidYMid meet"
     font-weight="$(canvas.fontweight)"
     font-family="$(canvas.fontfamily)"
     width="$(w)" height="$(h)"$(vb)>
    """
    bg = """<rect x="0" y="0" width="$(vbWf)" height="$(vbHf)"
     fill="#$(hex(canvas.bgcolor))" opacity="$(canvas.bgopacity)"/>
    """
    endsvg = "</svg>\n"
    return join([header, bg, canvas.bgelements..., canvas.elements..., endsvg], "")
end


"""
    drawsvg(mol::SimpleMolGraph) -> String

Generate molecular structure image as a SVG format string.

`width` and `height` specifies the size of the image (width and height
attribute of svg tag).
"""
function drawsvg(mol::SimpleMolGraph;
        bgcolor="#FFF", bgopacity=0.0,
        atomhighlight=eltype(mol)[], bondhighlight=Edge{eltype(mol)}[], highlightcolor="#FDD835",
        atomindex=false, indexcolor="#000", indexbgcolor="#F0F0FF",
        kwargs...)
    canvas = SvgCanvas(parse(RGB, bgcolor), bgopacity)
    draw2d!(canvas, mol; kwargs...)
    # highlight atoms if is_atom_visible=true or no incident edges
    # setdiff(Int[], []) -> Any[], setdiff(Int[], Int[]) -> Int[]  ???
    enodes = Set{eltype(mol)}(vcat([[src(e), dst(e)] for e in bondhighlight]...))
	nodes_to_show = collect(setdiff(
        atomhighlight, setdiff(enodes, findall(is_atom_visible(mol)))))
    sethighlight!(canvas, nodes_to_show, parse(RGB, highlightcolor))
    sethighlight!(canvas, bondhighlight, parse(RGB, highlightcolor))
    atomindex && drawatomindex!(canvas, is_atom_visible(mol), parse(RGB, indexcolor), parse(RGB, indexbgcolor))
    return tosvg(canvas; kwargs...)
end


"""
    html_fixed_size(mol::SimpleMolGraph, width, height) -> HTML{String}

Generate fixed-size HTML wrapper for the SVG element.

The width and height args can be numeric values (converted to `px`) or CSS strings like `100%`.
"""
function html_fixed_size(mol, width, height; kwargs...)
    wstr = width isa Real ? "$(convert(Int, round(width)))px" : width
    hstr = height isa Real ? "$(convert(Int, round(height)))px" : height
    svg = drawsvg(mol; kwargs...)
    return HTML("""<div style="width:$(wstr);height:$(hstr)">$(svg)</div>""")
end


"""
    html_grid(mols, cols, rowheight) -> HTML{String}

Generate grid layout HTML wrapper for the SVG elements.

`cols` - number of columns in the grid.
`rowheight` - numeric value (converted to `px`) or CSS string like `100%`.
"""
function html_grid(mols, cols, rowheight; kwargs...)
    htmls = String[]
    hstr = rowheight isa Real ? "$(convert(Int, round(rowheight)))px" : rowheight
    for row in Iterators.partition(mols, cols)
        push!(htmls, """<div style="display:grid; grid-template-columns:repeat($(cols), 1fr); grid-template-rows:1fr;">""")
        for m in row
            push!(htmls, drawsvg(m, height=hstr; kwargs...))
        end
        push!(htmls, "</div>")
    end
    return HTML(join(htmls, ""))
end


# Custom pretty printing for Plute notebook environment
Base.show(io::IO, m::MIME"text/html", mol::SimpleMolGraph
    ) = show(io, m, html_fixed_size(mol, "250px", "250px"))
Base.show(io::IO, m::MIME"text/html", mols::Vector{<:SimpleMolGraph}
    ) = show(io, m, html_grid(mols, 4, "200px"))

# Workaround until query visualizer implemented
Base.show(io::IO, m::MIME"text/html", mol::SimpleMolGraph{<:Integer,<:QueryTree,<:QueryTree}
    ) = print(io, "{$(nv(mol)), $(ne(mol))} simple molecular graph query $(typeof(mol))")
Base.show(io::IO, m::MIME"text/html", mols::Vector{<:SimpleMolGraph{<:Integer,<:QueryTree,<:QueryTree}}
    ) = print(io, "$(length(mols)) simple molecular graph queries $(eltype(mols))")


function initcanvas!(
        canvas::SvgCanvas, coords::AbstractArray{Float64}, boundary::Tuple)
    (top, left, width, height, unit) = boundary
    sf = canvas.scaleunit / unit
    pd = [canvas.paddingXf canvas.paddingYf] * canvas.scaleunit
    canvas.coords = (coords .- [left top]) .* [1 -1] .* sf .+ pd
    viewbox = ([width height] * sf) .+ (pd * 2)
    canvas.viewboxW = viewbox[1]
    canvas.viewboxH = viewbox[2]
    return
end


function atommarkupsvg(canvas::SvgCanvas, symbol, charge, implicith, direction)
    small = round(Int, canvas.fontsize * 0.7)
    return atommarkup(
        symbol, charge, implicith, direction,
        """<tspan baseline-shift="-25%" font-size="$(small)">""", "</tspan>",
        """<tspan baseline-shift="50%" font-size="$(small)">""", "</tspan>"
    )
end

atommarkupleft(canvas::SvgCanvas, symbol, charge, implicith) = atommarkupsvg(
    canvas, symbol, charge, implicith, :left)

atommarkupright(canvas::SvgCanvas, symbol, charge, implicith) = atommarkupsvg(
    canvas, symbol, charge, implicith, :right)


function drawtextsvg!(canvas::SvgCanvas, pos, text, color, anchor, xoffset)
    xy = svgcoords(pos + (xoffset, canvas.fontsize / 2))
    elem = """<text $(xy) font-size="$(canvas.fontsize)" fill="#$(hex(color))"$(anchor)>$(text)</text>"""
    push!(canvas.elements, elem)
    return
end

drawtextleft!(canvas::SvgCanvas, pos, text, color) = drawtextsvg!(
    canvas, pos, text, color, """ text-anchor="end" """, canvas.fontsize / 2)

drawtextcenter!(canvas::SvgCanvas, pos, text, color) = drawtextsvg!(
    canvas, pos, text, color, """ text-anchor="middle" """, 0)

drawtextright!(canvas::SvgCanvas, pos, text, color) = drawtextsvg!(
    canvas, pos, text, color, " ", canvas.fontsize / -2)


function drawtextannot!(canvas::SvgCanvas, pos, text, color, bgcolor)
    size = round(Int, canvas.fontsize * canvas.annotsizef)
    bxy = svgcoords(pos)
    txy = svgcoords(pos + (0, size))
    elem = """<g>
     <rect $(bxy) width="$(size)" height="$(size)" rx="$(size/2)" ry="$(size/2)" fill="#$(hex(bgcolor))" />
     <text $(txy) font-size="$(size)" fill="#$(hex(color))">$(text)</text>
    </g>
    """
    push!(canvas.elements, elem)
    return
end

function drawtexthighlight!(canvas::SvgCanvas, pos, color)
    size = round(Int, canvas.fontsize * canvas.hlsizef)
    xy = svgcoords(pos - (size / 2, size / 2))
    elem = """
     <rect $(xy) width="$(size)" height="$(size)" rx="$(size/2)" ry="$(size/2)" fill="#$(hex(color))" />
    """
    push!(canvas.bgelements, elem)
    return
end


function drawline!(canvas::SvgCanvas, seg, color; isdashed=false)
    option = isdashed ? """ stroke-dasharray="10,10" """ : " "
    coords = svgcoords(seg)
    elem = """<line $(coords) stroke="#$(hex(color))"$(option)/>"""
    push!(canvas.elements, elem)
    return
end

function drawline!(canvas::SvgCanvas, seg, ucolor, vcolor; isdashed=false)
    ucolor == vcolor && return drawline!(canvas, seg, ucolor, isdashed=isdashed)
    option = isdashed ? """ stroke-dasharray="10,10" """ : " "
    mid = midpoint(seg)
    coords1 = svgcoords(Segment(seg.u, mid))
    coords2 = svgcoords(Segment(mid, seg.v))
    elem = """<line $(coords1) stroke="#$(hex(ucolor))"$(option)/>
    <line $(coords2) stroke="#$(hex(vcolor))"$(option)/>
    """
    push!(canvas.elements, elem)
    return
end

drawdashedline!(canvas::SvgCanvas, seg, ucolor, vcolor) = drawline!(
    canvas, seg, ucolor, vcolor, isdashed=true)


function drawwedge!(canvas::SvgCanvas, seg, color)
    """ u ◀︎ v """
    d = distance(seg)
    scalef = Point2D(d, canvas.wedgewidthf / 2 * canvas.scaleunit)
    rotatef = unitvector(seg)
    translf = seg.u
    svgtf = svgtransform(transformmatrix(scalef, rotatef, translf))
    elem = """<polygon points="0,0 1,1 1,-1" fill="#$(hex(color))" transform="matrix($(svgtf))"/>
    """  # length: 1, width: 2
    push!(canvas.elements, elem)
    return
end

function drawwedge!(canvas::SvgCanvas, seg, ucolor, vcolor)
    """ u ◀︎ v """
    ucolor == vcolor && return drawwedge!(canvas, seg, ucolor)
    d = distance(seg)
    scalef = Point2D(d, canvas.wedgewidthf / 2 * canvas.scaleunit)
    rotatef = unitvector(seg)
    translf = seg.u
    svgtf = svgtransform(transformmatrix(scalef, rotatef, translf))
    elem = """<g stroke-width="0.3" transform="matrix($(svgtf))">
     <polygon points="0,0 0.5,-0.5 0.5,0.5" fill="#$(hex(ucolor))"/>
     <polygon points="0.5,-0.5 0.5,0.5 1,1 1,-1" fill="#$(hex(vcolor))"/>
     </g>
    """  # length: 1, width: 2
    push!(canvas.elements, elem)
    return
end


function drawdashedwedge!(canvas::SvgCanvas, seg, color)
    """ u ◁ v """
    d = distance(seg)
    scalef = Point2D(d / 7, canvas.wedgewidthf / 16 * canvas.scaleunit)
    rotatef = unitvector(seg)
    translf = seg.u
    svgtf = svgtransform(transformmatrix(scalef, rotatef, translf))
    elem = """<g stroke="#$(hex(color))" stroke-width="0.3" transform="matrix($(svgtf))">
     <line x1="0" y1="1" x2="0" y2="-1" />
     <line x1="1" y1="2" x2="1" y2="-2" />
     <line x1="2" y1="3" x2="2" y2="-3" />
     <line x1="3" y1="4" x2="3" y2="-4" />
     <line x1="4" y1="5" x2="4" y2="-5" />
     <line x1="5" y1="6" x2="5" y2="-6" />
     <line x1="6" y1="7" x2="6" y2="-7" />
     <line x1="7" y1="8" x2="7" y2="-8" />
    </g>
    """  # length: 7, width: 16
    push!(canvas.elements, elem)
    return
end

function drawdashedwedge!(canvas::SvgCanvas, seg, ucolor, vcolor)
    """ u ◁ v """
    ucolor == vcolor && return drawdashedwedge!(canvas, seg, ucolor)
    d = distance(seg)
    scalef = Point2D(d / 7, canvas.wedgewidthf / 16 * canvas.scaleunit)
    rotatef = unitvector(seg)
    translf = seg.u
    svgtf = svgtransform(transformmatrix(scalef, rotatef, translf))
    elem = """<g stroke-width="0.3" transform="matrix($(svgtf))">
     <line x1="0" y1="1" x2="0" y2="-1" stroke="#$(hex(ucolor))" />
     <line x1="1" y1="2" x2="1" y2="-2" stroke="#$(hex(ucolor))" />
     <line x1="2" y1="3" x2="2" y2="-3" stroke="#$(hex(ucolor))" />
     <line x1="3" y1="4" x2="3" y2="-4" stroke="#$(hex(ucolor))" />
     <line x1="4" y1="5" x2="4" y2="-5" stroke="#$(hex(vcolor))" />
     <line x1="5" y1="6" x2="5" y2="-6" stroke="#$(hex(vcolor))" />
     <line x1="6" y1="7" x2="6" y2="-7" stroke="#$(hex(vcolor))" />
     <line x1="7" y1="8" x2="7" y2="-8" stroke="#$(hex(vcolor))" />
    </g>
    """  # length: 7, width: 16
    push!(canvas.elements, elem)
    return
end


function drawwave!(canvas::SvgCanvas, seg, color)
    d = distance(seg)
    scalef = Point2D(d / 7, canvas.wedgewidthf / 2 * canvas.scaleunit)
    rotatef = unitvector(seg)
    translf = seg.u
    svgtf = svgtransform(transformmatrix(scalef, rotatef, translf))
    elem = """<polyline points="0,0 0.5,0 1,1 2,-1 3,1 4,-1 5,1 6,-1 6.5,0 7,0"
     stroke="#$(hex(color))" stroke-width="0.2" fill="none" transform="matrix($(svgtf))"/>
    """  # length: 7, width: 2
    push!(canvas.elements, elem)
    return
end

function drawwave!(canvas::SvgCanvas, seg, ucolor, vcolor)
    ucolor == vcolor && return drawwave!(canvas, seg, ucolor)
    d = distance(seg)
    scalef = Point2D(d / 7, canvas.wedgewidthf / 2 * canvas.scaleunit)
    rotatef = unitvector(seg)
    translf = seg.u
    svgtf = svgtransform(transformmatrix(scalef, rotatef, translf))
    elem = """<g stroke-width="0.2" fill="none" transform="matrix($(svgtf))">
     <polyline points="0,0 0.5,0 1,1 2,-1 3,1 3.5,0" stroke="#$(hex(ucolor))" />
     <polyline points="3.5,0 4,-1 5,1 6,-1 6.5,0 7,0" stroke="#$(hex(vcolor))" />
     </g>
    """  # length: 7, width: 2
    push!(canvas.elements, elem)
    return
end


function drawlinehighlight!(canvas::SvgCanvas, seg, color)
    cds = svgcoords(seg)
    w = canvas.linehlwidthf * canvas.scaleunit
    elem = """<line $(cds) stroke="#$(hex(color))" stroke-width="$(w)" stroke-linecap="round"/>
    """
    push!(canvas.bgelements, elem)
    return
end
