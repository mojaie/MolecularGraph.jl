#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    Color


struct Color
    r::Int
    g::Int
    b::Int
end


Formatting.format(expr::String, c::Color) = format(expr, c.r, c.g, c.b)
