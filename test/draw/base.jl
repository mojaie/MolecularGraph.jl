#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

using MolecularGraph:
    Color, atom_color, is_atom_visible,
    sdf_bond_style, double_bond_style, atomhtml

@testset "draw.base" begin
    @testset "atomnotation" begin
        ASSET_DIR = joinpath(dirname(@__FILE__), "..", "..", "assets", "test")
        mol = sdftomol(open(joinpath(ASSET_DIR, "demo.mol")))

        @test atom_color(mol)[6] == Color(0, 192, 0)
        @test atom_color(mol)[23] == Color(0, 192, 192)  # Default color
        @test is_atom_visible(mol)[1] === false
        @test is_atom_visible(mol)[8] === true
    end

    @testset "bondstyle" begin
        ASSET_DIR = joinpath(dirname(@__FILE__), "..", "..", "assets", "test")
        mol = sdftomol(open(joinpath(ASSET_DIR, "demo.mol")))
        sb = sdf_bond_style(mol)
        @test sb[7] === :up  # 3 ◀ 7
        @test sb[9] === :down  # 4 ◁ 8
        @test sb[10] === :unspecified # 5 ~ 11
        db = double_bond_style(mol.graph, bond_order(mol), coords2d(mol), sssr(mol))
        @test db[1] === :clockwise  # 2 ニ 1
        @test db[18] === :anticlockwise # 17 ニ 18
        @test db[12] === :none  # 10 = 13
    end

    @testset "atomhtml" begin
        @test atomhtml(:C, 0, 3, :right) == "CH<sub>3</sub>"
        @test atomhtml(:N, 1, 4, :left) == "<sup>+</sup>H<sub>4</sub>N"
    end
end
