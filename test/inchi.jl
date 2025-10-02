
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
    inchi3_no_stereo = inchi(mol_from_inchi, options = "SNon")
    # reconstruction from InChI now includes stereo information, but coordinates are
    # generated via coordgen!, which leads to different stereo perception in some cases,
    # e.g. definite parity instead of undefined
    @test_broken inchi3 == inchi1
    # without stereo, the test passes
    @test inchi3_no_stereo == inchi1_no_stereo
    # however, another round-trip via inchitomol should retain the stereo information
    mol2_from_inchi = inchitomol(inchi3)
    inchi4 = inchi(mol2_from_inchi)
    @test inchi3 == inchi4
    @test has_exact_match(mol_from_inchi, mol2_from_inchi)
    # Moreover, applying coordgen! to the original molecule also yields the same InChI
    coordgen!(mol)
    @test inchi(mol) == inchi3
end
