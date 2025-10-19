
@testset "json" begin

@testset "json.sdfile" begin
    atoms = [SDFAtom(),SDFAtom(),SDFAtom()]
    bonds = [SDFBond(),SDFBond()]
    mol = MolGraph(Edge.([(1, 2), (2, 3)]), atoms, bonds)
    mol2 = mol_from_json(JSON.json(mol))
    @test mol2 isa SDFMolGraph
    @test mol == mol2
    @test mol !== mol2

    atoms = [SMILESAtom(),SMILESAtom(),SMILESAtom()]
    bonds = [SMILESBond(),SMILESBond()]
    mol = MolGraph(Edge.([(1, 2), (2, 3)]), atoms, bonds)
    mol2 = mol_from_json(JSON.json(mol))
    @test mol2 isa SMILESMolGraph
    @test mol == mol2
    @test mol !== mol2

    assetdir = joinpath(dirname(@__FILE__), "..", "assets", "test")
    mol = sdftomol(joinpath(assetdir, "demo.mol"))
    mol2 = mol_from_json(JSON.json(mol))
    @test mol == mol2
    @test mol !== mol2
    mol = sdftomol(joinpath(assetdir, "aspirin_v3.mol"))
    mol2 = mol_from_json(JSON.json(mol))
    @test mol == mol2
    @test mol !== mol2

    # isolated vertices
    mol = sdftomol(joinpath(assetdir, "hydrate.mol"))
    mol2 = mol_from_json(JSON.json(mol))
    @test mol == mol2
    @test mol !== mol2

    # multibyte properties
    mol = sdftomol(joinpath(assetdir, "biotin_2d.sdf"))
    mol2 = mol_from_json(JSON.json(mol))
    @test mol == mol2
    @test mol !== mol2
    @test mol2["日本国内法規制情報"] == "該当なし"
end

@testset "json.smarts" begin
    atoms = [
        QueryAtom([(1, 2), (1, 3)], [qor(), qeq(:symbol, "N"), qeq(:symbol, "O")]),
        QueryAtom([(1, 2), (1, 3)], [qand(), qeq(:symbol, "C"), qeq(:isotope, "14")]),
        QueryAtom([(1, 2)], [qnot(), qtrue(:isaromatic)])
    ]
    bonds = [
        QueryBond(Tuple{Int,Int}[], [qeq(:order, "2")]),
        QueryBond([(1, 2)], [qnot(), qtrue(:isaromatic)])
    ]
    mol = QueryMolGraph(Edge.([(1, 2), (2, 3)]), atoms, bonds)
    mol2 = mol_from_json(JSON.json(mol))
    @test mol == mol2
    @test mol !== mol2

    mol = smartstomol(raw"[$([CX3]([#6])[#6]),$([CX3H][#6])]=[$([NX2][#6]),$([NX2H])]")
    mol2 = mol_from_json(JSON.json(mol))
    @test mol == mol2
    @test mol !== mol2

    mol = smartstomol("[O,S]=P([O,S])([O,S])[O,S]")
    mol2 = mol_from_json(JSON.json(mol))
    @test mol == mol2
    @test mol !== mol2
end

@testset "json.smiles" begin
    mol = smilestomol("OCCc1c(C)[n+](=cs1)Cc2cnc(C)nc(N)2")
    mol2 = mol_from_json(JSON.json(mol))
    @test mol == mol2
    @test mol !== mol2

    pinene = smilestomol("CC1([C@H]2CCC(=C)[C@@H]1C2)C")
    pinene2 = mol_from_json(JSON.json(pinene))
    @test pinene == pinene2
    @test pinene !== pinene2

    sildenafil = smilestomol("O=S(=O)(N1CCN(C)CC1)c4cc(c2[nH]c(=O)c3n(C)nc(CCC)c3n2)c(OCC)cc4")
    sildenafil2 = mol_from_json(JSON.json(sildenafil))
    @test sildenafil == sildenafil2
    @test sildenafil !== sildenafil2

    # isolated vertices
    hydrate = smilestomol("[Cu+2].[O-]S(=O)(=O)[O-].O.O.O.O.O")
    hydrate2 = mol_from_json(JSON.json(hydrate))
    @test hydrate == hydrate2
    @test hydrate !== hydrate2
end

@testset "json.rdkit" begin
    assetdir = joinpath(dirname(@__FILE__), "..", "assets", "test")
    mol = sdftomol(joinpath(assetdir, "demo.mol"))
    rmol = mol_from_json(to_rdkjson(mol))
    @test is_ring_aromatic(mol) == is_ring_aromatic(rmol)
    # @test smiles(mol) == smiles(rmol)
    mol2 = smilestomol("CN1CC[C@]23c4c5ccc(c4O[C@H]2[C@H](C=C[C@H]3[C@H]1C5)O)OC")
    rmol2 = mol_from_json(to_rdkjson(mol2))
    @test is_ring_aromatic(mol2) == is_ring_aromatic(rmol2)
    # @test smiles(mol2) == smiles(rmol2)

    # isolated vertices
    hydrate = smilestomol("[Cu+2].[O-]S(=O)(=O)[O-].O.O.O.O.O")
    hydrate2 = mol_from_json(to_rdkjson(hydrate))
    @test atom_counter(hydrate) == atom_counter(hydrate2)

    # 3D conformer
    asp = sdftomol(joinpath(assetdir, "aspirin_3d.sdf"))
    d = to_rdkdict(asp)
    @test length(d["molecules"][1]["conformers"][1]["coords"]) == 21

    # stereobond remapping
    nata = sdftomol(joinpath(assetdir, "nata.mol"))
    remove_all_hydrogens!(nata)
    d = to_rdkdict(nata)
    @test d["molecules"][1]["bonds"][32]["stereoAtoms"] == [24, 36]
end

end # json
