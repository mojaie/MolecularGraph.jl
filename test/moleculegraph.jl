@testset "moleculargraph" begin
    mol = MoleculeGraph()

    @test mol.atom[3].symbol == "C"
    @test mol.mw == 12.011
end
