
@testset "inchi" begin
    demomol = joinpath(dirname(@__FILE__), "..", "assets", "test", "demo.mol")
    txt = read(open(demomol), String)
    inchi1 = inchi(txt)
    ikey1 = inchikey(inchi1)
    demomol = joinpath(dirname(@__FILE__), "..", "assets", "test", "demo.mol")
    mol = sdftomol(demomol)
    inchi2 = inchi(mol)
    ikey2 = inchikey(inchi2)
    @test inchi1 == inchi2
    @test ikey1 == ikey2
end
