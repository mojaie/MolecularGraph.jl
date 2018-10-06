

@testset "sdfilereader" begin

"    1.1763    0.6815    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0"
"   -0.1041    2.3896    0.0000 Br  0  0  0  0  0  0  0  0  0  0  0  0"
"    1.6514    2.7627    0.0000 I   0  0  0  0  0  0  0  0  0  0  0  0"
"    2.5488    1.2083    0.0000 F   0  3  0  0  0  0  0  0  0  0  0  0"
"   -0.2917    0.6045    0.0000 N   0  5  0  0  0  0  0  0  0  0  0  0"
    atomblockpath = joinpath(
        dirname(@__FILE__), "..", "assets", "test", "atomblock1.txt")
    atomblock = readlines(atomblockpath)
    atoms = GraphMol.parseatom(atomblock)
    print(atoms)
    @test weight(atoms) == 12.011
end
