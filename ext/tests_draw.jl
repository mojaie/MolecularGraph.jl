#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

module MolecularGraphDrawTest

using MolecularGraph
using Test

import Cairo

@testset "draw2d" begin
    ASSETS_DIR = joinpath(dirname(@__FILE__), "..", "assets")
    OUT_DIR = joinpath(dirname(@__FILE__), "..", "_temp")
    mkpath(OUT_DIR)

    function draws(infile, outfile)
        m = sdftomol(open(joinpath(ASSETS_DIR, "test", infile)))
        io = open(joinpath(OUT_DIR, outfile), "w")
        write(io, drawsvg(m))
        close(io)
    end
    draws("demo.mol", "demo.svg")
    draws("biotin_2d.sdf", "biotin_2d.svg")
    draws("DB00115.mol", "DB00115.svg")
    draws("nata.mol", "nata.svg")
    draws("null.mol", "null.svg")

    function drawp(infile, outfile)
        m = sdftomol(open(joinpath(ASSETS_DIR, "test", infile)))
        io = open(joinpath(OUT_DIR, outfile), "w")
        drawpng(io, m, 400, 400)
        close(io)
    end
    drawp("demo.mol", "demo.png")
    drawp("biotin_2d.sdf", "biotin_2d.png")
    drawp("DB00115.mol", "DB00115.png")
    drawp("nata.mol", "nata.png")
    drawp("null.mol", "null.png")
end

end
