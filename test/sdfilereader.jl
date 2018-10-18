
@testset "sdfilereader" begin

@testset "parseatoms" begin
    atomblock = [
        "    1.1763    0.6815    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",
        "   -0.1041    2.3896    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0",
        "    1.6514    2.7627    0.0000 I   0  0  0  0  0  0  0  0  0  0  0  0",
        "    2.5488    1.2083    0.0000 F   0  3  0  0  0  0  0  0  0  0  0  0",
        "   -0.2917    0.6045    0.0000 N   0  5  0  0  0  0  0  0  0  0  0  0"
    ]
    atoms = parseatoms(atomblock)
    @test length(atoms) == 5
    @test atoms[1].index == 1
    @test atoms[2].symbol == "Br"
    @test atoms[3].coords == (1.6514f0, 2.7627f0, 0.0f0)
end

@testset "parsebonds" begin
    bondblock = [
        "  1  2  2  0  0  0  0",
        "  2  6  1  0  0  0  0",
        "  5  4  1  6  0  0  0"
    ]
    bonds = parsebonds(bondblock)
    @test length(bonds) == 3
    @test bonds[1].u == 1
    @test bonds[1].v == 2
    @test bonds[1].order == 2
    @test bonds[3].u == 4
    @test bonds[3].v == 5
    @test bonds[3].notation == 4
end

@testset "parsemol" begin
    demomol = joinpath(dirname(@__FILE__), "..", "assets", "test", "demo.mol")
    mol = parsemol(readlines(demomol))
    @test length(mol.graph.nodes) == 37
    @test length(mol.graph.edges) == 37
    @test length(mol.graph.adjacency) == 37
end

end # sdfilereader
