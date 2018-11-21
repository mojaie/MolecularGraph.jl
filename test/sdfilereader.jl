
@testset "sdfilereader" begin

@testset "parsesdfatom" begin
    full = parsesdfatom(
        "    1.1763    0.6815    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0"
    )
    @test full[1] == :C
    @test full[2] == 0
    @test full[3] == 1
    @test full[5] == [1.1763, 0.6815, 0.0]
    atomBr = parsesdfatom("   -0.1041    2.3896    0.0000 Br  0  0")
    @test atomBr[1] == :Br
    @test atomBr[5] == [-0.1041, 2.3896, 0.0]
    charged = parsesdfatom("    2.5488    1.2083    0.0000 O   0  3")
    @test charged[2] == 1
    radical = parsesdfatom("    1.6514    2.7627    0.0000 C   0  4")
    @test radical[2] == 0
    @test radical[3] == 2
end

@testset "parsesdfmolbond" begin
    bond1 = parsesdfbond("  1  2  2  0  0  0  0")
    @test bond1 == [1, 2, 2, 0]
    bond2 = parsesdfbond("  5  4  1  6  0  0  0")
    @test bond2 == [5, 4, 1, 6]
end

@testset "parsesdfmol" begin
    demomol = joinpath(dirname(@__FILE__), "..", "assets", "test", "demo.mol")
    mol = parsesdfmol(readlines(demomol))
    @test length(mol.graph.edges) == 37
    @test length(mol.graph.adjacency) == 37
    @test mol isa SDFile
end

end # sdfilereader
