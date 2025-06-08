
@testset "json" begin

@testset "sdfile" begin
    atoms = [SDFAtom(),SDFAtom(),SDFAtom()]
    bonds = [SDFBond(),SDFBond()]
    mol = MolGraph(Edge.([(1, 2), (2, 3)]), atoms, bonds)
    mol2 = MolGraph(to_json(mol))
    @test mol2 isa SDFMolGraph
    @test mol == mol2
    @test mol !== mol2

    atoms = [SMILESAtom(),SMILESAtom(),SMILESAtom()]
    bonds = [SMILESBond(),SMILESBond()]
    mol = MolGraph(Edge.([(1, 2), (2, 3)]), atoms, bonds)
    mol2 = MolGraph(to_json(mol))
    @test mol2 isa SMILESMolGraph
    @test mol == mol2
    @test mol !== mol2
end

@testset "serialization" begin
    mol = smilestomol("OCCc1c(C)[n+](=cs1)Cc2cnc(C)nc(N)2")
    mol2 = MolGraph(to_json(mol))
    @test mol == mol2
    @test mol !== mol2

    mol = smilestomol("C1=CC=C2C(=C1)C(=O)OC23C4=C(C=C(C=C4)O)OC5=C3C=CC(=C5)O")
    mol2 = MolGraph(to_json(mol))
    @test mol == mol2
    @test mol !== mol2

    sildenafil = smilestomol("O=S(=O)(N1CCN(C)CC1)c4cc(c2[nH]c(=O)c3n(C)nc(CCC)c3n2)c(OCC)cc4")
    sildenafil2 = MolGraph(to_json(sildenafil))
    @test sildenafil == sildenafil2
    @test sildenafil !== sildenafil2
end

@testset "smarts" begin
    mol = smartstomol(raw"[$([CX3]([#6])[#6]),$([CX3H][#6])]=[$([NX2][#6]),$([NX2H])]")
    mol2 = MolGraph(to_json(mol))
    @test mol == mol2
    @test mol !== mol2

    mol = smartstomol("[O,S]=P([O,S])([O,S])[O,S]")
    mol2 = MolGraph(to_json(mol))
    @test mol == mol2
    @test mol !== mol2
end

@testset "rdkit" begin
    assetdir = joinpath(dirname(@__FILE__), "..", "assets", "test")
    mol = sdftomol(joinpath(assetdir, "demo.mol"))
    rmol = MolGraph(to_rdkdict(mol))
    @test smiles(mol) == smiles(rmol)
    mol2 = smilestomol("CN1CC[C@]23c4c5ccc(c4O[C@H]2[C@H](C=C[C@H]3[C@H]1C5)O)OC")
    rmol2 = MolGraph(to_rdkdict(mol2))
    @test smiles(mol2) == smiles(rmol2)
end

end # json
