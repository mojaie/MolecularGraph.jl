
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
    @test full[:intrinsic_valence] == -1
    @test atomBr[:intrinsic_valence] == -1
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
    @test mols[2]["PubChem CID"] == "68827"

    sdf = joinpath(assetdir, "sdfile_error.sdf")
    mols = collect(sdfilereader(sdf))
    @test nv(mols[2]) == 0
    @test mols[2]["Name"] == "Artemisinin"

    sdf = joinpath(assetdir, "marvin_dative.sdf")
    mols = collect(sdfilereader(sdf))
    @test nv(mols[1]) == 0
    @test mols[1]["Name"] == "Auranofin"

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
    pcharges = split(mol["PUBCHEM_MMFF94_PARTIAL_CHARGES"], "\n")
    @test length(pcharges) == parse(Int, pcharges[1]) + 1
    for i = 2:length(pcharges)
        node, charge = split(pcharges[i])
        @test 1 <= parse(Int, node) <= nv(mol)
        @test isa(parse(Float32, charge), Float32)
    end
    coord_type = split(mol["PUBCHEM_COORDINATE_TYPE"])
    @test coord_type == ["2", "5", "10"]  # ensure the last attribute is read

end

end # sdfilereader

@testset "hydrogen flags of atom lines" begin
    flags_v2 = MolecularGraph.ctab_hydrogen_flags_v2
    flags_v3 = MolecularGraph.ctab_hydrogen_flags_v3
    atom(hhh, vvv) = "    0.0000    0.0000    0.0000 C   0  0  0" * lpad(hhh, 3) * "  0" * lpad(vvv, 3) * "  0  0  0  0  0  0"
    @test flags_v2(atom(0, 0)) == (false, 0)
    @test flags_v2(atom(1, 0)) == (true, 0)
    @test flags_v2(atom(0, 2)) == (false, 2)
    @test flags_v2(atom(0, 15)) == (false, -1)  # zero valence
    @test flags_v2("    1.6514    2.7627    0.0000 C   0  4") == (false, 0)  # short line
    @test flags_v3("M  V30 1 N 0 0 0 0") == (false, 0)
    @test flags_v3("M  V30 1 N 0 0 0 0 HCOUNT=-1") == (true, 0)
    @test flags_v3("M  V30 1 N 0 0 0 0 CHG=1 VAL=4") == (false, 4)
    @test flags_v3("M  V30 1 N 0 0 0 0 VAL=-1") == (false, -1)

    intrinsic_valence = MolecularGraph.intrinsic_valence_from_flags
    @test intrinsic_valence(false, 0) == -1  # automatic
    @test intrinsic_valence(false, 3) == 3  # the valence is the intrinsic valence
    @test intrinsic_valence(false, -1) == 0  # zero valence
    @test intrinsic_valence(true, 0) == 0  # hydrogen count is read as "no automatic hydrogens"
    @test intrinsic_valence(true, 3) == 3  # the valence wins
end

@testset "fixed number of hydrogens from SDFile" begin
    molblock(sym, hhh, vvv; chg="") = "\n  test\n\n  1  0  0  0  0  0  0  0  0  0999 V2000\n" *
        "    0.0000    0.0000    0.0000 $(rpad(sym, 3)) 0  0  0$(lpad(hhh, 3))  0$(lpad(vvv, 3))  0  0  0  0  0  0\n$(chg)M  END\n"
    rd(args...; kwargs...) = sdftomol(IOBuffer(molblock(args...; kwargs...)))
    @test implicit_hydrogens(rd("C", 0, 0)) == [4]  # unchanged
    @test implicit_hydrogens(rd("C", 0, 15)) == [0]  # charcoal
    @test implicit_hydrogens(rd("C", 1, 0)) == [0]  # hydrogen count flag
    @test implicit_hydrogens(rd("C", 0, 2)) == [2]
    @test rd("C", 0, 15)[1].intrinsic_valence == 0
    # hydroxide: charge from the property block must not reset the hydrogens
    @test implicit_hydrogens(rd("O", 0, 1; chg="M  CHG  1   1  -1\n")) == [1]
    @test implicit_hydrogens(rd("O", 0, 15; chg="M  CHG  1   1  -1\n")) == [0]
    # two atoms, bonded: total valence 3 of N with one single bond gives two hydrogens
    # (a valence below the apparent valence is like the apparent valence: no implicit hydrogens)
    mol = sdftomol(IOBuffer("\n  test\n\n  2  1  0  0  0  0  0  0  0  0999 V2000\n" *
        "    0.0000    0.0000    0.0000 N   0  0  0  0  0  3  0  0  0  0  0  0\n" *
        "    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" *
        "  1  2  1  0\nM  END\n"))
    @test [mol[i].intrinsic_valence for i in vertices(mol)] == [3, -1]
    @test implicit_hydrogens(mol) == [2, 3]
    low = sdftomol(IOBuffer("\n  test\n\n  3  2  0  0  0  0  0  0  0  0999 V2000\n" *
        "    0.0000    0.0000    0.0000 C   0  0  0  0  0  1  0  0  0  0  0  0\n" *
        "    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" *
        "    2.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n" *
        "  1  2  1  0\n  1  3  1  0\nM  END\n"))
    @test implicit_hydrogens(low) == [0, 3, 3]
    # removing and adding hydrogens keeps the number of hydrogens
    remove_hydrogens!(mol)
    @test implicit_hydrogens(mol) == [2, 3]
    add_hydrogens!(mol)
    @test nv(mol) == 7
    remove_hydrogens!(mol)
    @test implicit_hydrogens(mol) == [2, 3]
end

@testset "intrinsic valence above 4" begin
    # central atom with a valence column and `nb` single bonds to carbons
    function hyper(sym, vvv, nb)
        atoms = "    0.0000    0.0000    0.0000 $(rpad(sym, 3)) 0  0  0  0  0$(lpad(vvv, 3))  0  0  0  0  0  0\n"
        for k in 1:nb
            atoms *= "    $(k).0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0\n"
        end
        bonds = join(["  1$(lpad(k + 1, 3))  1  0\n" for k in 1:nb])
        sdftomol(IOBuffer("\n  test\n\n$(lpad(nb + 1, 3))$(lpad(nb, 3))  0  0  0  0  0  0  0  0999 V2000\n" * atoms * bonds * "M  END\n"))
    end
    sulfur = hyper("S", 6, 2)
    @test sulfur[1].intrinsic_valence == 6
    @test valence(sulfur)[1] == 6
    @test implicit_hydrogens(sulfur)[1] == 4
    @test lone_pair(sulfur)[1] == 0
    phosphorus = hyper("P", 5, 3)
    @test implicit_hydrogens(phosphorus)[1] == 2
    # without the valence column the element default applies
    @test implicit_hydrogens(hyper("S", 0, 2))[1] == 0
    @test implicit_hydrogens(hyper("N", 3, 1))[1] == 2
    # the valence column is written back and removing hydrogens keeps the valence
    @test split(printv2mol(sulfur), "\n")[5][49:51] == "  6"
    add_hydrogens!(sulfur)
    @test nv(sulfur) == 13  # 4 hydrogens on S, 3 on each of the carbons
    @test implicit_hydrogens(sulfur)[1] == 0
    remove_hydrogens!(sulfur)
    @test implicit_hydrogens(sulfur)[1] == 4
end
