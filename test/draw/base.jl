#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "draw.base" begin
    @testset "atomnotation" begin
        ASSET_DIR = joinpath(dirname(@__FILE__), "..", "..", "assets", "test")
        mol = sdftomol(open(joinpath(ASSET_DIR, "demo.mol")))

        @test atom_color(mol)[6] == RGB(0, 192/255, 0)
        @test atom_color(mol)[23] == RGB(0, 192/255, 192/255)  # Default color
        @test is_atom_visible(mol)[1] === false
        @test is_atom_visible(mol)[8] === true
    end

    @testset "bondstyle" begin
        ASSET_DIR = joinpath(dirname(@__FILE__), "..", "..", "assets", "test")
        mol = sdftomol(open(joinpath(ASSET_DIR, "demo.mol")))
        sb = draw2d_bond_style(mol)
        @test sb[7] === :up  # 3 ◀ 7
        @test sb[9] === :down  # 4 ◁ 8
        @test sb[10] === :unspecified # 5 ~ 11
        db = double_bond_style(mol.graph, bond_order(mol), coords2d(mol), sssr(mol))
        @test db[1] === :clockwise  # 2 ニ 1
        @test db[18] === :anticlockwise # 17 ニ 18
        @test db[12] === :none  # 10 = 13
    end

    @testset "atom_markup" begin
        @test atom_markup(:Br, 0, 0) == [[(:default, "Br")]]
        @test atom_markup(:N, 0, 4) == [[(:default, "N")], [(:default, "H"), (:sub, "4")]]
        @test atom_markup(:C, 13, 0) == [[(:sup, "13"), (:default, "C")]]
    end

    @testset "coords" begin
        coords = Coords2d([Point2d(x...) for x in [
            [3.6902, -1.0041],
            [3.6891, -1.8315],
            [4.4039, -2.2444],
            [5.1203, -1.831]
        ]])
        g = SimpleGraph([Edge(1, 2), Edge(2,3), Edge(3,4)])
        coords2 = normalize_coords(
            g, coords, zeros(Int, 4), fill(:right, 4), falses(4), 1.0, 1.0, 0.0, 0.0
        )
        @test all(isapprox.(extrema([p[1] for p in coords2[1]]), (0.0, 1.73034012866)))
        @test all(isapprox.(extrema([p[2] for p in coords2[1]]), (0.0, 1.49953945052)))
        @test isapprox.(coords2[2], 1.73034012866)
        @test isapprox.(coords2[3], 1.49953945052)
    end
end
