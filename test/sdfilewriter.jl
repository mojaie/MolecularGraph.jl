
@testset "sdfilewriter" begin

@testset "molblock" begin
    demomol = joinpath(dirname(@__FILE__), "..", "assets", "test", "demo.mol")
    mol = sdftomol(demomol)
    dump = printv2mol(mol)
    mol2 = sdftomol(IOBuffer(dump))
    @test nv(mol2) == 37

    stereohyd = joinpath(dirname(@__FILE__), "..", "assets", "test", "DB00115.mol")
    mol = sdftomol(stereohyd)
    dump = printv2mol(mol; givebackhydrogen=false)  # TODO: VirtualAtom would be better
    mol2 = sdftomol(IOBuffer(dump))
    @test nv(mol2) == 93

    mol3d = joinpath(dirname(@__FILE__), "..", "assets", "test", "aspirin_3d.sdf")
    mol = sdftomol(mol3d)
    dump = printv2mol(mol)
    mol2 = sdftomol(IOBuffer(dump))
    @test length(coords3d(mol2)) == 21

    smol = smilestomol("CCC1CC(C=O)CCC1N")
    smol2 = sdftomol(IOBuffer(printv2mol(smol)))
    @test nv(smol2) == 11
    @test isempty(smol2[:metadata])
end

end #sdfilewriter
