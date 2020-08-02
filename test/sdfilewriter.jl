
@testset "sdfilewriter" begin

@testset "molblock" begin
    demomol = joinpath(dirname(@__FILE__), "..", "assets", "test", "demo.mol")
    mol = sdftomol(demomol)
    mol2 = sdftomol(split(printv2mol(mol), "\n"))
    @test nodecount(mol2) == 37
    
    smol = smilestomol("CCC1CC(C=O)CCC1N")
    smol2 = sdftomol(split(printv2mol(smol), "\n"))
    @test nodecount(smol2) == 11
end

end #sdfilewriter
