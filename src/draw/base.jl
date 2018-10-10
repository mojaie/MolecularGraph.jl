#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#


export
    Canvas

abstract type Canvas end

CANVAS = Dict(:svg => SvgCanvas)

function draw(format::Symbol, mol::MolecularGraph)
    canvas = CANVAS[format]()
    draw!(canvas, mol)
    show(canvas)
