#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    SvgCanvas,
    drawline!,
    drawwedge!

mutable struct SvgCanvas <: Canvas
    fontweight::String
    fontfamily::String
    bgcolor::Color
    opacity::UInt8
    header::Vector{String}
    components::Vector{String}

    function SvgCanvas(symbol::AbstractString)
        canvas = new()
        canvas.fontweight = "normal"
        canvas.fontfamily = "Helvetica"
        canvas.bgcolor = Color("#ffffff")
        canvas.opacity = 0

        canvas.mbwidthf = 0.15
        canvas.triminnerf = 0.2
        canvas.trimoverlapf = 0.3

        canvas.header = [
            '<svg xmlns="http://www.w3.org/2000/svg"',
            ' xmlns:xlink="http://www.w3.org/1999/xlink"',
            ' version="1.2" baseProfile="tiny"',
            ' text-rendering="geometricPrecision"'
        ]
        canvas.components = []
        canvas
    end
end


function coordsconv(canvas::Canvas, p::Point2D)
    canvas.size / 2 + vinv(p) * canvas.scalefactor
end


function drawline!(canvas::Canvas, u::Point2D, v::Point2D, color::Color;
                   vcolor=nothing, dashed=false)
    option = dashed ? """ stroke-dasharray="10,10" """ : " "
    if vcolor !== nothing && color != vcolor
        mid = u + v / 2
        push!(canvas.components,
              """<line x1="$(u.x)" y1="$(u.y)" x2="$(mid.x)" y2="$(mid.y)"
               stroke="$(color)"$(option)/>\n""")
        push!(canvas.components,
              """<line x1="$(mid.x)" y1="$(mid.y)" x2="$(v.x)" y2="$(v.y)"
               stroke="$(vcolor)"$(option)/>\n""")
    else
        push!(canvas.components,
              """<line x1="$(u.x)" y1="$(u.y)" x2="$(v.x)" y2="$(v.y)"
               stroke="$()"$(option)/>\n""")
    end
end


function drawwedge!(canvas::Canvas, head::Point2D, tail::Point2D, color::Color)
    width = displaymlb * 0.2
    t1 = parallelmove(tail, head, -2pi, width / 2)[1]
    t2 = parallelmove(tail, head, 2pi, width / 2)[1]
    vrts = join(format(v, 2) for v in [head, t1, t2], "")
    push!(canvas.components,
          """<polygon points="$(vrts)" fill="$(color)" />\n""")
