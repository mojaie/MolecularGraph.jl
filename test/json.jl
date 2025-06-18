
# @testset "json" begin

@testset "json.sdfile" begin
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

    assetdir = joinpath(dirname(@__FILE__), "..", "assets", "test")
    mol = sdftomol(joinpath(assetdir, "demo.mol"))
    mol2 = MolGraph(to_json(mol))
    @test mol == mol2
    @test mol !== mol2
    mol = sdftomol(joinpath(assetdir, "aspirin_v3.mol"))
    mol2 = MolGraph(to_json(mol))
    @test mol == mol2
    @test mol !== mol2
end

@testset "json.smarts" begin
    atoms = [
        QueryAtom([(1, 2), (1, 3)], [qor(), qeq(:symbol, "N"), qeq(:symbol, "O")]),
        QueryAtom([(1, 2), (1, 3)], [qand(), qeq(:symbol, "C"), qeq(:mass, "14")]),
        QueryAtom([(1, 2)], [qnot(), qtrue(:isaromatic)])
    ]
    bonds = [
        QueryBond(Tuple{Int,Int}[], [qeq(:order, "2")]),
        QueryBond([(1, 2)], [qnot(), qtrue(:isaromatic)])
    ]
    mol = QueryMolGraph(Edge.([(1, 2), (2, 3)]), atoms, bonds)
    mol2 = QueryMolGraph(to_json(mol))
    @test mol == mol2
    @test mol !== mol2

    mol = smartstomol(raw"[$([CX3]([#6])[#6]),$([CX3H][#6])]=[$([NX2][#6]),$([NX2H])]")
    mol2 = QueryMolGraph(to_json(mol))
    @test mol == mol2
    @test mol !== mol2

    mol = smartstomol("[O,S]=P([O,S])([O,S])[O,S]")
    mol2 = QueryMolGraph(to_json(mol))
    @test mol == mol2
    @test mol !== mol2
end

@testset "json.smiles" begin
    mol = smilestomol("OCCc1c(C)[n+](=cs1)Cc2cnc(C)nc(N)2")
    mol2 = MolGraph(to_json(mol))
    @test mol == mol2
    @test mol !== mol2

    pinene = smilestomol("CC1([C@H]2CCC(=C)[C@@H]1C2)C")
    pinene2 = MolGraph(to_json(pinene))
    @test pinene == pinene2
    @test pinene !== pinene2

    sildenafil = smilestomol("O=S(=O)(N1CCN(C)CC1)c4cc(c2[nH]c(=O)c3n(C)nc(CCC)c3n2)c(OCC)cc4")
    sildenafil2 = MolGraph(to_json(sildenafil))
    @test sildenafil == sildenafil2
    @test sildenafil !== sildenafil2
end

@testset "json.rdkit" begin
    assetdir = joinpath(dirname(@__FILE__), "..", "assets", "test")
    mol = sdftomol(joinpath(assetdir, "demo.mol"))
    rmol = MolGraph(to_rdkdict(mol))
    @test smiles(mol) == smiles(rmol)
    mol2 = smilestomol("CN1CC[C@]23c4c5ccc(c4O[C@H]2[C@H](C=C[C@H]3[C@H]1C5)O)OC")
    rmol2 = MolGraph(to_rdkdict(mol2))
    @test smiles(mol2) == smiles(rmol2)
end

# end # json
