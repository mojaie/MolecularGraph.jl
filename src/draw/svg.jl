#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#


mutable struct SvgCanvas <: Canvas
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

    viewboxW::Float64
    viewboxH::Float64
    bgelements::Vector{String}
    elements::Vector{String}
    coords::Vector{Point2d}
    optpos::Vector{Point2d}
    optdeg::Vector{Float64}

    function SvgCanvas(bgcolor::RGB, bgopacity::Float64)
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
        canvas.hlsizef = 1.2
        canvas.paddingXf = 1.0
        canvas.paddingYf = 1.0

        # Appearance
        canvas.fontweight = "normal"
        canvas.fontfamily = "Helvetica"
        canvas.fontsize = 14.0
        canvas.fonttagmap = Dict(
            :sub => ("""<tspan baseline-shift="-25%" font-size="70%">""", "</tspan>"),
            :sup => ("""<tspan baseline-shift="50%" font-size="70%">""", "</tspan>")
        )
        canvas.bgcolor = bgcolor
        canvas.bgopacity = bgopacity

        # Canvas state
        canvas.viewboxW = 1.0
        canvas.viewboxH = 1.0
        canvas.bgelements = []
        canvas.elements = []
        canvas.optpos = []
        canvas.optdeg = []

        return canvas
    end
end


svgcolor(c::RGB) = @sprintf "rgb(%d, %d, %d)" c.r*255 c.g*255 c.b*255
svgcoords(p::Point) = @sprintf "x=\"%.2f\" y=\"%.2f\"" p[1] p[2]
svgcirclecoords(p::Point) = @sprintf "cx=\"%.2f\" cy=\"%.2f\"" p[1] p[2]
svgcoords(s::Line
    ) = @sprintf "x1=\"%.2f\" y1=\"%.2f\" x2=\"%.2f\" y2=\"%.2f\"" s[1][1] s[1][2] s[2][1] s[2][2]
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
     fill="$(svgcolor(canvas.bgcolor))" opacity="$(canvas.bgopacity)"/>
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
        bgcolor="#FFF", bgopacity=1.0,
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
function html_fixed_size(mol::SimpleMolGraph, width, height; kwargs...)
    wstr = width isa Real ? "$(convert(Int, round(width)))px" : width
    hstr = height isa Real ? "$(convert(Int, round(height)))px" : height
    svg = drawsvg(mol; kwargs...)
    return HTML("""<div style="width:$(wstr);height:$(hstr)">$(svg)</div>""")
end


"""
    html_grid(mols, cols::Int, rowheight) -> HTML{String}

Generate grid layout HTML wrapper for the SVG elements.

`cols` - number of columns in the grid.
`rowheight` - numeric value (converted to `px`) or CSS string like `100%`.
"""
function html_grid(mols, cols::Int, rowheight; kwargs...)
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
Base.show(io::IO, m::MIME"text/html", mol::QueryMolGraph
    ) = print(io, "{$(nv(mol)), $(ne(mol))} simple molecular graph query $(typeof(mol))")
Base.show(io::IO, m::MIME"text/html", mols::Vector{QueryMolGraph}
    ) = print(io, "$(length(mols)) simple molecular graph queries $(eltype(mols))")


function initcanvas!(
        canvas::SvgCanvas, coords::Vector{Point2d}, boundary::Tuple)
    (top, left, width, height, unit) = boundary
    sf = canvas.scaleunit / unit
    pd = [canvas.paddingXf canvas.paddingYf] * canvas.scaleunit
    conv = p -> (p - Point2d(left, top)) * Point2d(1, -1) * sf + Point2d(pd...)
    canvas.coords = conv.(coords)
    viewbox = ([width height] * sf) .+ (pd * 2)
    canvas.viewboxW = viewbox[1]
    canvas.viewboxH = viewbox[2]
    return
end


const SVG_ATOM_TEXT_ANCHOR = Dict(
        :left => """ text-anchor="end" """, :center => """ text-anchor="middle" """, :right => "")
const SVG_ATOM_TEXT_XOFFSET = Dict(:left => 1, :center => 0, :right => -1)


function drawtext!(canvas::SvgCanvas, pos, text, color, align)
    anchor = SVG_ATOM_TEXT_ANCHOR[align]
    xoffset = SVG_ATOM_TEXT_XOFFSET[align]
    # empirical offset factor (0.4, 0.35)
    xy = svgcoords(pos + Point2d(canvas.fontsize * 0.4 * xoffset, canvas.fontsize * 0.35))
    elem = """<text $(xy) font-size="$(canvas.fontsize)" fill="$(svgcolor(color))"$(anchor)>$(text)</text>"""
    push!(canvas.elements, elem)
    return
end


function drawtextannot!(canvas::SvgCanvas, pos, text, color, bgcolor)
    size = round(Int, canvas.fontsize * canvas.annotsizef)
    bxy = svgcoords(pos)
    txy = svgcoords(pos + Point2d(0, size))
    elem = """<g>
     <rect $(bxy) width="$(size)" height="$(size)" rx="$(size/2)" ry="$(size/2)" fill="$(svgcolor(bgcolor))" />
     <text $(txy) font-size="$(size)" fill="$(svgcolor(color))">$(text)</text>
    </g>
    """
    push!(canvas.elements, elem)
    return
end

function drawtexthighlight!(canvas::SvgCanvas, pos, color)
    size = round(Int, canvas.fontsize * canvas.hlsizef)
    xy = svgcoords(pos - Point2d(size / 2, size / 2))
    elem = """
     <rect $(xy) width="$(size)" height="$(size)" rx="$(size/2)" ry="$(size/2)" fill="$(svgcolor(color))" />
    """
    push!(canvas.bgelements, elem)
    return
end


function drawline!(canvas::SvgCanvas, seg, color; isdashed=false)
    option = isdashed ? """ stroke-dasharray="10,10" """ : " "
    coords = svgcoords(seg)
    elem = """<line $(coords) stroke="$(svgcolor(color))"$(option)/>"""
    push!(canvas.elements, elem)
    return
end

function drawline!(canvas::SvgCanvas, seg, ucolor, vcolor; isdashed=false)
    ucolor == vcolor && return drawline!(canvas, seg, ucolor, isdashed=isdashed)
    option = isdashed ? """ stroke-dasharray="10,10" """ : " "
    mid = (seg[1] + seg[2]) / 2
    coords1 = svgcoords(Line(seg[1], mid))
    coords2 = svgcoords(Line(mid, seg[2]))
    elem = """<line $(coords1) stroke="$(svgcolor(ucolor))"$(option)/>
    <line $(coords2) stroke="$(svgcolor(vcolor))"$(option)/>
    """
    push!(canvas.elements, elem)
    return
end

drawdashedline!(canvas::SvgCanvas, seg, ucolor, vcolor) = drawline!(
    canvas, seg, ucolor, vcolor, isdashed=true)


function drawwedge!(canvas::SvgCanvas, seg, color)
    """ u ◀︎ v """
    svgtf = svgtransform(transformmatrix(seg, 1.0, canvas.wedgewidth))
    elem = """<polygon points="0,0 1,1 1,-1" fill="$(svgcolor(color))" transform="matrix($(svgtf))"/>
    """  # length: 1, width: 2
    push!(canvas.elements, elem)
    return
end

function drawwedge!(canvas::SvgCanvas, seg, ucolor, vcolor)
    """ u ◀︎ v """
    ucolor == vcolor && return drawwedge!(canvas, seg, ucolor)
    svgtf = svgtransform(transformmatrix(seg, 1.0, canvas.wedgewidth))
    elem = """<g stroke-width="0.3" transform="matrix($(svgtf))">
     <polygon points="0,0 0.5,-0.5 0.5,0.5" fill="$(svgcolor(ucolor))"/>
     <polygon points="0.5,-0.5 0.5,0.5 1,1 1,-1" fill="$(svgcolor(vcolor))"/>
     </g>
    """  # length: 1, width: 2
    push!(canvas.elements, elem)
    return
end


function drawdashedwedge!(canvas::SvgCanvas, seg, color)
    """ u ◁ v """
    svgtf = svgtransform(transformmatrix(seg, 1 / 7, canvas.wedgewidth / 8))
    elem = """<g stroke="$(svgcolor(color))" stroke-width="0.3" transform="matrix($(svgtf))">
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
    svgtf = svgtransform(transformmatrix(seg, 1 / 7, canvas.wedgewidth / 8))
    elem = """<g stroke-width="0.3" transform="matrix($(svgtf))">
     <line x1="0" y1="1" x2="0" y2="-1" stroke="$(svgcolor(ucolor))" />
     <line x1="1" y1="2" x2="1" y2="-2" stroke="$(svgcolor(ucolor))" />
     <line x1="2" y1="3" x2="2" y2="-3" stroke="$(svgcolor(ucolor))" />
     <line x1="3" y1="4" x2="3" y2="-4" stroke="$(svgcolor(ucolor))" />
     <line x1="4" y1="5" x2="4" y2="-5" stroke="$(svgcolor(vcolor))" />
     <line x1="5" y1="6" x2="5" y2="-6" stroke="$(svgcolor(vcolor))" />
     <line x1="6" y1="7" x2="6" y2="-7" stroke="$(svgcolor(vcolor))" />
     <line x1="7" y1="8" x2="7" y2="-8" stroke="$(svgcolor(vcolor))" />
    </g>
    """  # length: 7, width: 16
    push!(canvas.elements, elem)
    return
end


function drawwave!(canvas::SvgCanvas, seg, color)
    svgtf = svgtransform(transformmatrix(seg, 1 / 7, canvas.wedgewidth))
    elem = """<polyline points="0,0 0.5,0 1,1 2,-1 3,1 4,-1 5,1 6,-1 6.5,0 7,0"
     stroke="$(svgcolor(color))" stroke-width="0.2" fill="none" transform="matrix($(svgtf))"/>
    """  # length: 7, width: 2
    push!(canvas.elements, elem)
    return
end

function drawwave!(canvas::SvgCanvas, seg, ucolor, vcolor)
    ucolor == vcolor && return drawwave!(canvas, seg, ucolor)
    svgtf = svgtransform(transformmatrix(seg, 1 / 7, canvas.wedgewidth))
    elem = """<g stroke-width="0.2" fill="none" transform="matrix($(svgtf))">
     <polyline points="0,0 0.5,0 1,1 2,-1 3,1 3.5,0" stroke="$(svgcolor(ucolor))" />
     <polyline points="3.5,0 4,-1 5,1 6,-1 6.5,0 7,0" stroke="$(svgcolor(vcolor))" />
     </g>
    """  # length: 7, width: 2
    push!(canvas.elements, elem)
    return
end


function drawlinehighlight!(canvas::SvgCanvas, seg, color)
    cds = svgcoords(seg)
    w = canvas.linehlwidth
    elem = """<line $(cds) stroke="$(svgcolor(color))" stroke-width="$(w)" stroke-linecap="round"/>
    """
    push!(canvas.bgelements, elem)
    return
end
