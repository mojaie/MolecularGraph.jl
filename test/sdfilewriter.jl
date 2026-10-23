
@testset "sdfilewriter" begin

@testset "molblock" begin
    demomol = joinpath(dirname(@__FILE__), "..", "assets", "test", "demo.mol")
    mol = sdftomol(demomol)
    dump = printv2mol(mol)
    mol2 = sdftomol(IOBuffer(dump))
    @test nv(mol2) == 37

    # stereospecific hydrogens are given back without regenerating 2D coordinates
    biotin = sdftomol(joinpath(dirname(@__FILE__), "..", "assets", "test", "biotin_2d.sdf"))
    ikey = inchikey(biotin)
    remove_all_hydrogens!(biotin)
    crds = collect(coords2d(biotin))
    gb = sdftomol(IOBuffer(printv2mol(biotin)))
    @test nv(gb) == nv(biotin) + 2  # ring fusion
    @test all(isapprox(collect(coords2d(gb))[i], crds[i]; atol=1e-3) for i in 1:nv(biotin))
    @test inchikey(gb) == ikey
    # 3D coordinates are kept
    aspirin = sdftomol(joinpath(dirname(@__FILE__), "..", "assets", "test", "aspirin_3d.sdf"))
    @test occursin("3D", split(printv2mol(aspirin), "\n")[2])
    # header line 2: program name (8 characters) with the package version
    ver = MolecularGraph.VERSION
    program = "MGjlv" * string(ver.major)[end] * lpad(string(ver.minor), 2, '0')
    @test split(printv2mol(aspirin), "\n")[2][3:10] == program

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

    # fixed number of hydrogens is written to the valence column and survives a round trip
    hmol = sdftomol(IOBuffer("\n  test\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n" *
        "    0.0000    0.0000    0.0000 N   0  0  0  0  0  3  0  0  0  0  0  0\n" *
        "    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" *
        "    2.0000    0.0000    0.0000 O   0  0  0  0  0  2  0  0  0  0  0  0\n" *
        "  1  2  1  0\n  2  3  1  0\nM  END\n"))
    @test [hmol[i].implicit_hs for i in vertices(hmol)] == [2, -1, 1]
    hdump = printv2mol(hmol)
    @test [l[49:51] for l in split(hdump, "\n")[5:7]] == ["  3", "  0", "  2"]
    @test [l[43:45] for l in split(hdump, "\n")[5:7]] == ["  0", "  0", "  0"]
    hmol2 = sdftomol(IOBuffer(hdump))
    @test [hmol2[i].implicit_hs for i in vertices(hmol2)] == [2, -1, 1]
    # zero valence
    cmol = sdftomol(IOBuffer("\n  test\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n    0.0000    0.0000    0.0000 C   0  0  0  0  0 15  0  0  0  0  0  0\nM  END\n"))
    cdump = printv2mol(cmol)
    @test split(cdump, "\n")[5][49:51] == " 15"
    @test implicit_hydrogens(sdftomol(IOBuffer(cdump))) == [0]
    # automatic hydrogens are not marked
    @test all(l[49:51] == "  0" for l in split(printv2mol(smilestomol("CCO")), "\n")[5:7])

    # SMILES bracket atoms have a fixed number of hydrogens, which survives a round trip
    for smiles in ("[O-]", "[C]", "C[NH3+]", "[CH2]C")
        bmol = smilestomol(smiles)
        remove_hydrogens!(bmol)
        bmol2 = sdftomol(IOBuffer(printv2mol(bmol)))
        @test total_hydrogens(bmol2) == total_hydrogens(bmol)
    end
    @test split(printv2mol(smilestomol("[O-]")), "\n")[5][49:51] == " 15"  # as RDKit writes it
    @test split(printv2mol(smilestomol("[OH-]")), "\n")[5][49:51] == "  1"
    @test split(printv2mol(smilestomol("[C]")), "\n")[5][49:51] == " 15"
    # stereo hydrogens given back to atoms with a fixed number of hydrogens are not counted twice
    fused = smilestomol("[CH]1(F)CC[C@H]2CCCC[C@@H]12")
    remove_hydrogens!(fused)
    fused2 = sdftomol(IOBuffer(printv2mol(fused)))
    @test sum(total_hydrogens(fused2)) == sum(total_hydrogens(fused))

    sdf = joinpath(dirname(@__FILE__), "..", "assets", "test", "drugbank1-100.sdf")
    io = IOBuffer()
    sdfilewriter(io, sdfilereader(sdf), givebackhydrogen=false)
    seekstart(io)
    @test sum(nv(mol) for mol in sdfilereader(io)) == sum(nv(mol) for mol in sdfilereader(sdf))
end

end #sdfilewriter
