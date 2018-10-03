@testset "atom" begin
    atom = Atom("C")

    @test atom.symbol == "C"
    @test weight(atom) == 12.011
end
