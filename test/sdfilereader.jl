
@testset "sdfilereader" begin

@testset "ctab_atom_v2" begin
    full = ctab_atom_v2(SDFAtom,
        "    1.1763    0.6815    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0"
    )
    @test full[:symbol] == :C
    @test full[:charge] == 0
    @test full[:multiplicity] == 1
    @test full[:coords] == [1.1763, 0.6815, 0.0]
    atomBr = ctab_atom_v2(SDFAtom, "   -0.1041    2.3896    0.0000 Br  0  0")
    @test atomBr[:symbol] == :Br
    @test atomBr[:coords] == [-0.1041, 2.3896, 0.0]
    charged = ctab_atom_v2(SDFAtom, "    2.5488    1.2083    0.0000 O   0  3")
    @test charged[:charge] == 1
    radical = ctab_atom_v2(SDFAtom, "    1.6514    2.7627    0.0000 C   0  4")
    @test radical[:charge] == 0
    @test radical[:multiplicity] == 2
end

@testset "ctab_bond_v2" begin
    edge, prop = ctab_bond_v2(Int, SDFBond, "  1  2  2  0  0  0  0")
    @test edge == Edge(1, 2)
    @test prop[:order] == 2
    @test prop[:notation] == 0
    @test prop[:isordered]
    edge, prop = ctab_bond_v2(Int, SDFBond, "  5  4  1  6  0  0  0")
    @test edge == Edge(4, 5)
    @test prop[:order] == 1
    @test prop[:notation] == 6
    @test !prop[:isordered]
end

@testset "sdftomol" begin
    assetdir = joinpath(dirname(@__FILE__), "..", "assets", "test")

    null = joinpath(assetdir, "null.mol")
    mol = sdftomol(null)
    @test ne(mol) == 0
    @test nv(mol) == 0
    @test mol isa SDFMolGraph

    demo = joinpath(assetdir, "demo.mol")
    mol = sdftomol(demo)
    @test ne(mol) == 37
    @test nv(mol) == 37

    sdf = joinpath(assetdir, "sdfile_test.sdf")
    mols = collect(sdfilereader(sdf))
    @test length(mols) == 3
    @test get_prop(mols[2], "PubChem CID") == "68827"

    sdf = joinpath(assetdir, "sdfile_error.sdf")
    mols = collect(sdfilereader(sdf))
    @test nv(mols[2]) == 0
    @test get_prop(mols[2], "Name") == "Artemisinin"

    aspirin_v3 = joinpath(assetdir, "aspirin_v3.mol")
    mol = sdftomol(aspirin_v3)
    @test ne(mol) == 13
    @test nv(mol) == 13

    kinase = joinpath(assetdir, "RHEA17826.rxn")
    rxn = rxntoreaction(kinase)
    @test length(rxn.reactants) == 2
    @test length(rxn.products) == 3

    kinase_rd = joinpath(assetdir, "RHEA17826.rd")
    rxn = collect(rdfilereader(kinase_rd))[1]
    @test length(rxn.reactants) == 2
    @test length(rxn.products) == 3

    kinase_v3 = joinpath(assetdir, "RHEA17826_v3.rxn")
    rxn = rxntoreaction(kinase_v3)
    @test length(rxn.reactants) == 2
    @test length(rxn.products) == 3

    # attributes
    aspirin_3d = joinpath(assetdir, "aspirin_3d.sdf")
    mol = sdftomol(aspirin_3d)
    pcharges = split(get_prop(mol, "PUBCHEM_MMFF94_PARTIAL_CHARGES"), "\n")
    @test length(pcharges) == parse(Int, pcharges[1]) + 1
    for i = 2:length(pcharges)
        node, charge = split(pcharges[i])
        @test 1 <= parse(Int, node) <= nv(mol)
        @test isa(parse(Float32, charge), Float32)
    end
    coord_type = split(get_prop(mol, "PUBCHEM_COORDINATE_TYPE"))
    @test coord_type == ["2", "5", "10"]  # ensure the last attribute is read

end

end # sdfilereader
