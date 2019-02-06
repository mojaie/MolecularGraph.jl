
@testset "draw.base" begin
    @testset "atomnotation" begin
        ASSET_DIR = joinpath(dirname(@__FILE__), "..", "..", "assets", "test")
        mol = sdftomol(open(joinpath(ASSET_DIR, "demo.mol")))
        atomnotation2d!(mol)
        @test mol[:Color2D][6] == Color(0, 192, 0)
        @test mol[:Color2D][23] == Color(0, 192, 192)  # Default color
        @test mol[:Visible2D][1] === false
        @test mol[:Visible2D][8] === true
    end

    @testset "bondnotation" begin
        ASSET_DIR = joinpath(dirname(@__FILE__), "..", "..", "assets", "test")
        mol = sdftomol(open(joinpath(ASSET_DIR, "demo.mol")))
        # TODO: Requires Cartesian2D
        matrix = zeros(Float64, nodecount(mol), 2)
        for (i, node) in nodesiter(mol)
            matrix[i, :] = node.coords[1:2]
        end
        mol.coords[:Cartesian2D] = cartesian2d(matrix)
        bondnotation2d!(mol)
        @test mol[:BondNotation][1] == 0  # 2 ニ 1
        @test mol[:BondNotation][21] == 1 # 17 ニ 18
        @test mol[:BondNotation][8] == 2  # 10 = 13
        @test mol[:BondNotation][25] == 0 # 27 ≡ 28
        @test mol[:BondNotation][15] == 1 # 3 ◀ 7
        @test mol[:BondNotation][5] == 6  # 4 ◁ 8
        @test mol[:BondNotation][28] == 4 # 5 ~ 11
    end

    @testset "atomhtml" begin
        @test atomhtml(:C, 0, 3, :right) == "CH<sub>3</sub>"
        @test atomhtml(:N, 1, 4, :left) == "<sup>+</sup>H<sub>4</sub>N"
    end
end
