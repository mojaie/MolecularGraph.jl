@testset "moleculargraph" begin
    mol = MolecularGraph()

    @test mol.atom[3].symbol == "C"
    @test mol.mw == 12.011
end
