#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

# CairoExt (requires Cairo). See also ext/tests_draw.jl for visual inspection of the outputs.

@testset "ext.cairo" begin
    TEST_DIR = joinpath(dirname(@__FILE__), "..", "..", "assets", "test")
    PNG_SIGNATURE = UInt8[0x89, 0x50, 0x4e, 0x47, 0x0d, 0x0a, 0x1a, 0x0a]
    function png(file; kwargs...)
        io = IOBuffer()
        drawpng(io, sdftomol(joinpath(TEST_DIR, file)), 400, 400; kwargs...)
        return take!(io)
    end

    @testset "drawpng" begin
        for file in ("demo.mol", "biotin_2d.sdf", "DB00115.mol", "nata.mol", "null.mol")
            @test png(file)[1:8] == PNG_SIGNATURE
        end
    end

    @testset "options" begin
        data = png("aspirin_3d.sdf";
            atomhighlight=[1, 2], bondhighlight=[Edge(1, 5)], atomindex=true, bgcolor="#F0F0F0")
        @test data[1:8] == PNG_SIGNATURE
        @test data != png("aspirin_3d.sdf")
        @test_throws ArgumentError png("aspirin_3d.sdf"; bondhighlight=[Edge(1, 2)])  # no bond
    end
end
