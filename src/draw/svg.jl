#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#


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


function tosvg(canvas::SvgCanvas; width="100%", height="100%", kwargs...)
    vbWf = @sprintf "%.2f" canvas.viewboxW
    vbHf = @sprintf "%.2f" canvas.viewboxH
    bgc = svgcolor(canvas.bgcolor)
    header = """<svg xmlns="http://www.w3.org/2000/svg"
     xmlns:xlink="http://www.w3.org/1999/xlink"
     version="1.1" baseProfile="tiny"
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
    return join([header, bg, canvas.bgelements..., canvas.elements..., endsvg], "")
end


"""
    drawsvg(mol::SimpleMolGraph) -> String

Generate molecular structure image as a SVG format string.

`width` and `height` specifies the size of the image (width and height
attribute of svg tag).
"""
function drawsvg(mol::SimpleMolGraph;
        atomhighlight=eltype(mol)[], bondhighlight=Edge{eltype(mol)}[], highlightcolor=Color(253, 216, 53),
        atomindex=false, indexcolor=Color(0, 0, 0), indexbgcolor=Color(240, 240, 255),
        kwargs...)
    canvas = SvgCanvas()
    draw2d!(canvas, mol; kwargs...)
    # highlight atoms if is_atom_visible=true or no incident edges
    # setdiff(Int[], []) -> Any[], setdiff(Int[], Int[]) -> Int[]  ???
    enodes = Set{eltype(mol)}(vcat([[src(e), dst(e)] for e in bondhighlight]...))
	nodes_to_show = collect(setdiff(
        atomhighlight, setdiff(enodes, findall(is_atom_visible(mol)))))
    sethighlight!(canvas, nodes_to_show, highlightcolor)
    sethighlight!(canvas, bondhighlight, highlightcolor)
    atomindex && drawatomindex!(canvas, is_atom_visible(mol), indexcolor, indexbgcolor)
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
    c = svgcolor(color)
    elem = """<text $(xy) font-size="$(canvas.fontsize)" fill="$(c)"$(anchor)>$(text)</text>"""
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
    c = svgcolor(color)
    bc = svgcolor(bgcolor)
    elem = """<g>
     <rect $(bxy) width="$(size)" height="$(size)" rx="$(size/2)" ry="$(size/2)" fill="$(bc)" />
     <text $(txy) font-size="$(size)" fill="$(c)">$(text)</text>
    </g>
    """
    push!(canvas.elements, elem)
    return
end

function drawtexthighlight!(canvas::SvgCanvas, pos, color)
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
    d = distance(seg)
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
    d = distance(seg)
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
    d = distance(seg)
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
    d = distance(seg)
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
    d = distance(seg)
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
    d = distance(seg)
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


function drawlinehighlight!(canvas::SvgCanvas, seg, color)
    cds = svgcoords(seg)
    c = svgcolor(color)
    elem = """<line $(cds) stroke="$(c)" stroke-width="10" stroke-linecap="round"/>
    """
    push!(canvas.bgelements, elem)
    return
end
