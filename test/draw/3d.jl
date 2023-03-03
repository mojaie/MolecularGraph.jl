using MakieCore

@testset "draw3d" begin
    mol = sdftomol(joinpath(@__DIR__, "..", "..", "assets", "test", "aspirin_3d.sdf"))
    @test length(atom_radius(mol)) == nv(mol)
end
