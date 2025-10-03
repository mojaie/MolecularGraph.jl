
@testset "model.atom" begin

@testset "atom" begin
    @test atom_number(:Fe) == 26
    @test atom_symbol(118) === :Og

    sdfa = SDFAtom(;symbol=:Fe, charge=2, coords=[1.0, 2.0, 0.0])
    @test sdfa[:symbol] === :Fe
    @test sdfa[:charge] == 2
    @test sdfa[:multiplicity] == 1
    @test sdfa[:coords] == [1.0, 2.0, 0.0]
    sdfa2 = SDFAtom(;symbol=:Fe, charge=2, coords=[1.0, 2.0, 0.0])
    @test sdfa == sdfa2
    @test hash(sdfa) == hash(sdfa2)
    @test_throws ErrorException SDFAtom(;symbol=:Ga, isotope=270)

    smia = SMILESAtom(isotope=13, stereo=:clockwise)
    @test smia[:symbol] === :C
    @test smia[:charge] == 0
    @test smia[:multiplicity] == 1
    @test smia[:isotope] == 13
    @test !smia[:isaromatic]
    @test smia[:stereo] === :clockwise
    smia2 = SMILESAtom(isotope=13, stereo=:clockwise)
    @test smia == smia2
    @test hash(smia) == hash(smia2)
    @test_throws ErrorException SMILESAtom(;symbol=:Ga, isotope=270)

end

end # model.atom