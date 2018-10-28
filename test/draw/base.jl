
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

    @testset "terminal_double_bond" begin
        vec = zeros(Int, 4)
        adj = [
            Dict(2=>1), Dict(1=>1, 3=>2), Dict(2=>2, 4=>3),
            Dict(3=>3, 5=>4), Dict(4=>4)
        ]
        numb = [1, 2, 2, 2, 1]
        valence = [2, 3, 3, 3, 1]
        terminal_double_bond!(vec, adj, numb, valence)
        @test vec == [2, 0, 0, 0]
    end

    @testset "draw2d" begin
        ASSET_DIR = joinpath(dirname(@__FILE__), "..", "..", "assets", "test")
        demomol = loadsdfmol(open(joinpath(ASSET_DIR, "demo.mol")))
        draw2d_annot!(demomol)
        @test count(demomol.v[:AtomVisible]) == 15
        @test count(demomol.v[:BondNotation] .== 2) == 2
        @test count(demomol.v[:BondNotation] .== 2) == 2
    end
end
