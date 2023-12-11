
@testset "sdfilewriter" begin

@testset "molblock" begin
    demomol = joinpath(dirname(@__FILE__), "..", "assets", "test", "demo.mol")
    mol = sdftomol(demomol)
    dump = printv2mol(mol)
    mol2 = sdftomol(IOBuffer(dump))
    @test nv(mol2) == 37

    smol = smilestomol("CCC1CC(C=O)CCC1N")
    smol2 = sdftomol(IOBuffer(printv2mol(smol)))
    @test nv(smol2) == 11
    @test isempty(metadata(smol2))
end

end #sdfilewriter
