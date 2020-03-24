
@testset "atom" begin
    sdfa = SDFileAtom(:Fe, 2, 1, nothing, [1.0, 2.0, 0.0], :anticlockwise)
    @test sdfa.symbol === :Fe
    @test sdfa.charge == 2
    @test sdfa.coords == [1.0, 2.0, 0.0]
    @test sdfa.stereo === :anticlockwise
    @test atomnumber(sdfa) == 26
    @test atomweight(sdfa) == 55.845

    smia = SmilesAtom(:C, 0, 1, 13.0, false, :clockwise)
    @test !smia.isaromatic
    @test smia.stereo === :clockwise
    @test atomname(smia) == "Carbon"
    @test atomweight(smia) == 13.0

    @test atomnumber(:Ra) == 88
    @test atomname(:Ti) == "Titanium"
end
