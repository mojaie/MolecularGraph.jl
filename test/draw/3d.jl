using AbstractPlotting

@testset "draw3d" begin
    mol = sdftomol(joinpath(@__DIR__, "..", "..", "assets", "test", "aspirin_3d.sdf"))
    # These are just tests to see if the call succeeds; they do not test the output
    @test isa(spacefilling(mol), Scene)
    @test isa(spacefilling(mol; H=false), Scene)
    @test isa(spacefilling(mol; tform=x->2x), Scene)

    @test isa(ballstick(mol), Scene)
    @test isa(ballstick(mol; H=false), Scene)
    @test isa(ballstick(mol; tform=x->2x), Scene)
    @test isa(ballstick(mol; tform=x->2x, markersize=2), Scene)
end
