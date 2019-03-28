
@testset "draw.base" begin
    @testset "atomnotation" begin
        ASSET_DIR = joinpath(dirname(@__FILE__), "..", "..", "assets", "test")
        mol = sdftomol(open(joinpath(ASSET_DIR, "demo.mol")))

        @test mol[:atomcolor][6] == Color(0, 192, 0)
        @test mol[:atomcolor][23] == Color(0, 192, 192)  # Default color
        @test mol[:isatomvisible][1] === false
        @test mol[:isatomvisible][8] === true
    end

    @testset "bondnotation" begin
        ASSET_DIR = joinpath(dirname(@__FILE__), "..", "..", "assets", "test")
        mol = sdftomol(open(joinpath(ASSET_DIR, "demo.mol")))
        @test mol[:bondnotation][1] == 0  # 2 ニ 1
        @test mol[:bondnotation][21] == 1 # 17 ニ 18
        @test mol[:bondnotation][8] == 2  # 10 = 13
        @test mol[:bondnotation][25] == 0 # 27 ≡ 28
        @test mol[:bondnotation][15] == 1 # 3 ◀ 7
        @test mol[:bondnotation][5] == 6  # 4 ◁ 8
        @test mol[:bondnotation][28] == 4 # 5 ~ 11
    end

    @testset "atomhtml" begin
        @test atomhtml(:C, 0, 3, :right) == "CH<sub>3</sub>"
        @test atomhtml(:N, 1, 4, :left) == "<sup>+</sup>H<sub>4</sub>N"
    end
end
