
@testset "atom" begin
    sdfa = sdfatom(:Fe, 2, 1, nothing, SVector(1.0, 2.0, 0.0))
    @test sdfa.symbol == :Fe
    @test sdfa.charge == 2
    @test sdfa.sdf_coords == [1.0, 2.0, 0.0]
    @test atomnumber(sdfa) == 26
    @test atomweight(sdfa) == 55.845

    smia = smilesatom(:C, 0, 1, 13.0, false, 2)
    @test !smia.smiles_aromatic
    @test smia.smiles_stereo == 2
    @test atomname(smia) == "Carbon"
    @test atomweight(smia) == 13.0

    @test atomnumber(:Ra) == 88
    @test atomname(:Ti) == "Titanium"
end
