
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

    # specified valence (vvv) is written and read back
    carbene = smilestomol("Cl[C]Cl")
    dump = printv2mol(carbene)
    @test split(dump, "\n")[6][49:51] == "  2"
    carbene2 = sdftomol(IOBuffer(dump))
    @test atom_valence(carbene2) == [-1, 2, -1]
    @test molecular_formula(carbene2) == "CCl2"
    atomc = sdftomol(IOBuffer(printv2mol(smilestomol("[C]"))))
    @test atom_valence(atomc) == [0]
    # the valence is kept if the atom is overwritten by the property block
    chg = smilestomol("C")
    chg[1] = SMILESAtom(symbol=:C, charge=-1, multiplicity=2, valence=2)
    chg2 = sdftomol(IOBuffer(printv2mol(chg)))
    @test atom_valence(chg2) == [2]
    @test atom_charge(chg2) == [-1]
    @test multiplicity(chg2) == [2]

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
