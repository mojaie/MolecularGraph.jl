@testset "atom" begin
    atom = GraphMol.Atom("C")

    @test atom.symbol == "C"
    @test !atom.visible
end
