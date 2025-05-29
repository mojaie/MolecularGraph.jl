
@testset "json" begin

@testset "serialization" begin
    assetdir = joinpath(dirname(@__FILE__), "..", "assets", "test")
    mol = sdftomol(joinpath(assetdir, "demo.mol"))
    rmol = rdktomol(to_rdkdict(mol))
    @test moltosmiles(mol) == moltosmiles(rmol)
    mol2 = smilestomol("CN1CC[C@]23c4c5ccc(c4O[C@H]2[C@H](C=C[C@H]3[C@H]1C5)O)OC")
    rmol2 = rdktomol(to_rdkdict(mol2))
    @test moltosmiles(mol2) == moltosmiles(rmol2)
end

end # json
