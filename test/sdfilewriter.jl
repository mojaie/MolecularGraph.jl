
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
    @test split(dump, "\n")[2][21:22] == "3D"

    smol = smilestomol("CCC1CC(C=O)CCC1N")
    sdump = printv2mol(smol)
    @test split(sdump, "\n")[2][21:22] == "2D"
    smol2 = sdftomol(IOBuffer(sdump))
    @test nv(smol2) == 11
    @test isempty(smol2[:metadata])

    sdf = joinpath(dirname(@__FILE__), "..", "assets", "test", "drugbank1-100.sdf")
    io = IOBuffer()
    sdfilewriter(io, sdfilereader(sdf), givebackhydrogen=false)
    seekstart(io)
    @test sum(nv(mol) for mol in sdfilereader(io)) == sum(nv(mol) for mol in sdfilereader(sdf))
end

end #sdfilewriter
