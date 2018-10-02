@testset "atom" begin
    atom = Atom("C")

    @test atom.symbol == "C"
    @test atom.std_weight == 12.011
end
