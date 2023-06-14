
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
    mol_from_inchi = inchitomol(inchi1)
    inchi1_no_stereo = inchi(txt, options = "SNon")
    inchi3 = inchi(mol_from_inchi)
    # lbinchi version 1.05 doesn't support chiral sdf output, therefore test without stereo information
    # as soon as version 1.06 is in place, the lines below need to be adapted
    @test_broken inchi3 == inchi1
    @test inchi3 == inchi1_no_stereo
end
