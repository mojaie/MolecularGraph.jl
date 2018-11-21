
@testset "atom" begin
    sdfa = SDFileAtom(:Fe, 2, 1, nothing, SVector(1.0, 2.0, 0.0))
    @test sdfa.symbol == :Fe
    @test sdfa.charge == 2
    @test sdfa.coords == [1.0, 2.0, 0.0]
    @test atomnumber(sdfa) == 26
    @test atomweight(sdfa) == 55.845

    smia = SmilesAtom(:C, 0, 1, 13.0, false, 2)
    @test !smia.isaromatic
    @test smia.stereo == 2
    @test atomname(smia) == "Carbon"
    @test atomweight(smia) == 13.0

    @test atomnumber(:Ra) == 88
    @test atomname(:Ti) == "Titanium"
end
