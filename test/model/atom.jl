
@testset "atom" begin
    atom = Atom("C")
    @test atom.symbol == "C"
    @test !atom.visible
    atom.Hcount = 3
    @test html(atom, :right) == "CH<sub>3</sub>"
    @test html(atom, :left) == "H<sub>3</sub>C"
end
