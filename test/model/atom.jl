
@testset "atom" begin
    atom = Atom("C")

    @test atom.symbol == "C"
    @test !atom.visible
end
