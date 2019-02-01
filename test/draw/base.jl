
@testset "draw.base" begin
    @testset "atomcolor" begin
        @test atomcolor(DRAW_SETTING)(:N) == Color(0, 0, 255)
        @test atomcolor(DRAW_SETTING)(:Fe) == Color(0, 192, 192)
        customsetting = copy(DRAW_SETTING)
        customsetting[:atomcolor][:S] = Color(128, 192, 0)
        customsetting[:atomcolor][:Mg] = Color(0, 128, 128)
        @test atomcolor(customsetting)(:S) == Color(128, 192, 0)
        @test atomcolor(customsetting)(:Mg) == Color(0, 128, 128)
    end

    @testset "atomvisible" begin
        @test atomvisible(false)(:C, 0)
        @test !atomvisible(false)(:C, 1)
        @test !atomvisible(false)(:C, 2)
        @test atomvisible(true)(:C, 0)
        @test atomvisible(true)(:C, 1)
        @test !atomvisible(true)(:C, 2)
    end

    @testset "atomhtml" begin
        @test atomhtml(:C, 0, 3, :right) == "CH<sub>3</sub>"
        @test atomhtml(:N, 1, 4, :left) == "<sup>+</sup>H<sub>4</sub>N"
    end

    @testset "draw2d" begin
        ASSET_DIR = joinpath(dirname(@__FILE__), "..", "..", "assets", "test")
        demomol = sdftomol(open(joinpath(ASSET_DIR, "demo.mol")))
        draw2d_annot!(demomol)
        @test count(demomol[:AtomVisible]) == 15
        @test count(demomol[:BondNotation] .== 2) == 2
        @test count(demomol[:BondNotation] .== 2) == 2
    end
end
