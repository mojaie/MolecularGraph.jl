
@testset "model.atom" begin

@testset "atom" begin
    @test atomnumber(:Fe) == 26
    @test atomsymbol(118) === :Og

    sdfa = SDFAtom(:Fe, 2, 1, nothing, [1.0, 2.0, 0.0])
    @test sdfa[:symbol] === :Fe
    @test sdfa[:charge] == 2
    @test sdfa[:multiplicity] == 1
    @test sdfa[:coords] == [1.0, 2.0, 0.0]
    sdfa2 = SDFAtom(["Fe", 2, 1, nothing, [1.0, 2.0, 0.0]])
    @test hash(sdfa) == hash(sdfa2)
    @test sdfa == sdfa2

    smia = SMILESAtom(:C, 0, 1, 13.0, false, :clockwise)
    @test smia[:symbol] === :C
    @test smia[:charge] == 0
    @test smia[:multiplicity] == 1
    @test smia[:mass] == 13
    @test !smia[:isaromatic]
    @test smia[:stereo] === :clockwise
    smia2 = SMILESAtom(["C", 0, 1, 13.0, false, :clockwise])
    @test hash(smia) == hash(smia2)
    @test smia == smia2

end

end # model.atom