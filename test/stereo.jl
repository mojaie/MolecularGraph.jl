#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "stereo" begin

@testset "isclockwise" begin
    @test isclockwise((1, 3, 4, true), 1, 4, 5)  # 3-4 4-5 (same)
    @test isclockwise((1, 3, 4, false), 1, 3, 5)  # 4-5 (reverse)
    @test !isclockwise((1, 3, 6, true), 3, 6, 7)  # 1-3 1-6 1-7 (reverse)
    @test !isclockwise((2, 5, 6, false), 2, 6, 7)  # 5-6 5-7 (same)
end

@testset "property" begin
    prop = MolProperty{Int}()
    prop.stereocenter[4] = (2, 6, 7, true)
    prop.stereocenter[5] = (4, 8, 9, false)
    prop.stereocenter[10] = (8, 3, 11, false)
    remap!(
        Val(:stereocenter), prop, [11, 2, 3, 4, 10, 6, 7, 8, 9],
        Edge.([(2, 4), (4, 5)])  # Not used
    )
    @test length(prop.stereocenter) == 2
    @test prop.stereocenter[4] == (2, 6, 7, true)
    @test prop.stereocenter[5] == (8, 3, 1, false)
    dump = to_dict(Val(:stereocenter), Val(:default), prop)
    @test reconstruct(
        Val(:stereocenter), MolProperty{Int}, dump) == prop.stereocenter

    prop = MolProperty{Int}()
    prop.stereobond[Edge(2 => 3)] = (1, 5, true)
    prop.stereobond[Edge(7 => 10)] = (8, 11, false)
    remap!(
        Val(:stereobond), prop, [1, 11, 3, 4, 10, 6, 7, 8, 9],
        Edge.([(1, 2), (2, 3)])  # Not used
    )
    @test length(prop.stereobond) == 1
    @test prop.stereobond[Edge(5 => 7)] == (8, 2, false)
    dump = to_dict(Val(:stereobond), Val(:default), prop)
    @test reconstruct(
        Val(:stereobond), MolProperty{Int}, dump) == prop.stereobond
end

@testset "stereo_hydrogen" begin
    # safe_remove_hydrogen! is already applied
    LAla1 = smilestomol("N[C@@H](C)C(=O)O")
    @test LAla1[2].stereo === :clockwise
    @test LAla1[:stereocenter][2] == (1, 4, 5, true)  # No.3 -> :H

    LAla2 = smilestomol("N[C@](C)([H])C(=O)O")
    @test LAla2[2].stereo === :anticlockwise
    @test LAla2[:stereocenter][2] == (1, 3, 5, true)  # No.4 -> :H

    LAla3 = smilestomol("[H][C@@](C(=O)O)(C)N")
    @test LAla3[2].stereo === :clockwise
    @test LAla3[:stereocenter][2] == (3, 6, 7, false)  # No.1 -> :H

    LAla4 = smilestomol("OC(=O)[C@]([H])(C)N")
    @test LAla4[4].stereo === :anticlockwise
    @test LAla4[:stereocenter][4] == (2, 6, 7, false)  # No.5 -> :H

    atomoxetine = smilestomol("O(c1ccccc1C)[C@H](c2ccccc2)CCNC")
    @test atomoxetine[:stereocenter][9] == (1, 11, 17, false)
    rem_vertex!(atomoxetine, 10)  # :H can be safely removed
    @test  atomoxetine[:stereocenter][9] == (1, 11, 17, false)

    ldopa = smilestomol("c1c(O)c(o)ccc1C[C@@H](N)C(=O)O")
    @test ldopa[:stereocenter][10] == (9, 12, 13, true)
    subg, vmap = induced_subgraph(ldopa, collect(9:15))
    # revmap 9 => 1, 10 => 2, 11 => 3, 12 => 4, 13 => 5, 14 => 6, 15 => 7
    @test subg[:stereocenter][2] == (1, 4, 5, true)
    rem_vertices!(ldopa, collect(1:8))
    # revmap 9 => 7, 10 => 6, 11 => 5, 12 => 4, 13 => 3, 14 => 2, 15 => 1
    @test ldopa[:stereocenter][6] == (7, 4, 3, true)
end


@testset "stereocenter_from_smiles" begin
    # ring bond
    asc = smilestomol("C([C@@H]([C@@H]1C(=C(C(=O)O1)O)O)O)O")
    @test asc[:smarts_lexical_succ][2] == [3, 4, 13]
    @test asc[:smarts_lexical_succ][4] == [5, 10, 6]
    @test asc[:stereocenter][2] == (1, 4, 13, true)
    @test asc[:stereocenter][4] == (2, 6, 10, false)
    subg, vmap = induced_subgraph(asc, [2, 4, 5, 6, 10])
    @test subg[:smarts_lexical_succ][2] == [3, 5, 4]
    @test subg[:stereocenter][2] == (1, 4, 5, false)
    codeine = smilestomol("CN1CC[C@]23c4c5ccc(c4O[C@H]2[C@H](C=C[C@H]3[C@H]1C5)O)OC")
    @test codeine[:smarts_lexical_succ][5] == [13, 19, 6]
    @test codeine[:smarts_lexical_succ][19] == [20, 5, 21]
    @test codeine[:stereocenter][5] == (4, 13, 19, false)
    @test codeine[:stereocenter][19] == (5, 18, 21, true)
end


@testset "stereocenter_from_sdf2d" begin
    # degree=4
    # default_logger = global_logger(ConsoleLogger(stdout, Logging.Debug))
    atoms = [
        SDFAtom(;coords=[0.0, 0.0, 0.0]),
        SDFAtom(;coords=[-1.41, 0.5, 0.0]),
        SDFAtom(;coords=[1.41, 0.5, 0.0]),
        SDFAtom(;coords=[-0.5, -1.41, 0.0]),
        SDFAtom(;coords=[0.5, -1.41, 0.0])
    ]
    bonds = [SDFBond() for _ in 1:4]
    edges = Edge.([(1,2), (1,3), (1,4), (1,5)])
    uns = MolGraph(edges, atoms, bonds, on_init=sdf_on_init!)
    @test isempty(uns[:stereocenter])

    bonds = [
        SDFBond(;notation=0),
        SDFBond(;notation=0),
        SDFBond(;notation=1),
        SDFBond(;notation=6)
    ]
    mol1 = MolGraph(edges, atoms, bonds, on_init=sdf_on_init!)
    @test mol1[:stereocenter][1] == (2, 3, 5, false)

    bonds = [
        SDFBond(;notation=1),
        SDFBond(;notation=0),
        SDFBond(;notation=6),
        SDFBond(;notation=0)
    ]
    mol2 = MolGraph(edges, atoms, bonds, on_init=sdf_on_init!)
    @test mol2[:stereocenter][1] == (2, 3, 5, true)

    bonds = [
        SDFBond(;notation=1),
        SDFBond(;notation=0),
        SDFBond(;notation=0),
        SDFBond(;notation=1)
    ]
    mol3 = MolGraph(edges, atoms, bonds, on_init=sdf_on_init!)
    @test mol3[:stereocenter][1] == (2, 3, 5, true)

    bonds = [
        SDFBond(;notation=6),
        SDFBond(;notation=0),
        SDFBond(;notation=0),
        SDFBond(;notation=6)
    ]
    mol4 = MolGraph(edges, atoms, bonds, on_init=sdf_on_init!)
    @test mol4[:stereocenter][1] == (2, 3, 5, false)

    bonds = [
        SDFBond(;notation=0),
        SDFBond(;notation=1),
        SDFBond(;notation=0),
        SDFBond(;notation=0)
    ]
    mol5 = MolGraph(edges, atoms, bonds, on_init=sdf_on_init!)
    @test mol5[:stereocenter][1] == (2, 3, 5, false)

    bonds = [
        SDFBond(;notation=0),
        SDFBond(;notation=0),
        SDFBond(;notation=6),
        SDFBond(;notation=0)
    ]
    mol6 = MolGraph(edges, atoms, bonds, on_init=sdf_on_init!)
    @test mol6[:stereocenter][1] == (2, 3, 5, true)

    bonds = [
        SDFBond(;notation=6),
        SDFBond(;notation=6),
        SDFBond(;notation=6),
        SDFBond(;notation=0)
    ]
    mol7 = MolGraph(edges, atoms, bonds, on_init=sdf_on_init!)
    @test mol7[:stereocenter][1] == (2, 3, 5, true)

    bonds = [
        SDFBond(;notation=1),
        SDFBond(;notation=0),
        SDFBond(;notation=1),
        SDFBond(;notation=1)
    ]
    mol8 = MolGraph(edges, atoms, bonds, on_init=sdf_on_init!)
    @test mol8[:stereocenter][1] == (2, 3, 5, true)

    bonds = [
        SDFBond(;notation=1),
        SDFBond(;notation=0),
        SDFBond(;notation=6),
        SDFBond(;notation=1)
    ]
    mol9 = MolGraph(edges, atoms, bonds, on_init=sdf_on_init!)
    @test mol9[:stereocenter][1] == (2, 3, 5, true)

    bonds = [
        SDFBond(;notation=6),
        SDFBond(;notation=1),
        SDFBond(;notation=1),
        SDFBond(;notation=6)
    ]
    mol10 = MolGraph(edges, atoms, bonds, on_init=sdf_on_init!)
    @test mol10[:stereocenter][1] == (2, 3, 5, false)

    bonds = [
        SDFBond(;notation=1),
        SDFBond(;notation=0),
        SDFBond(;notation=0),
        SDFBond(;notation=6)
    ]
    wrong1 = MolGraph(edges, atoms, bonds, on_init=sdf_on_init!)
    @test isempty(wrong1[:stereocenter])
    @test haskey(wrong1[:logs], "warning_stereocenter_ignored")
    # serialization check
    @test nv(MolGraph(to_json(wrong1))) == 5

    bonds = [
        SDFBond(;notation=1),
        SDFBond(;notation=1),
        SDFBond(;notation=0),
        SDFBond(;notation=0)
    ]
    wrong2 = MolGraph(edges, atoms, bonds, on_init=sdf_on_init!)
    @test isempty(wrong2[:stereocenter])
    @test haskey(wrong2[:logs], "warning_stereocenter_ignored")

    bonds = [
        SDFBond(;notation=1),
        SDFBond(;notation=1),
        SDFBond(;notation=1),
        SDFBond(;notation=1)
    ]
    wrong3 = MolGraph(edges, atoms, bonds, on_init=sdf_on_init!)
    @test isempty(wrong3[:stereocenter])
    @test haskey(wrong3[:logs], "warning_stereocenter_ignored")

    bonds = [
        SDFBond(;notation=6),
        SDFBond(;notation=6),
        SDFBond(;notation=6),
        SDFBond(;notation=6)
    ]
    wrong4 = MolGraph(edges, atoms, bonds, on_init=sdf_on_init!)
    @test isempty(wrong4[:stereocenter])
    @test haskey(wrong4[:logs], "warning_stereocenter_ignored")

    bonds = [
        SDFBond(;notation=1),
        SDFBond(;notation=6),
        SDFBond(;notation=1),
        SDFBond(;notation=6)
    ]
    wrong5 = MolGraph(edges, atoms, bonds, on_init=sdf_on_init!)
    @test isempty(wrong5[:stereocenter])
    @test haskey(wrong5[:logs], "warning_stereocenter_ignored")

    # degree=3
    atoms = [
        SDFAtom(;coords=[0.0, 0.0, 0.0]),
        SDFAtom(;coords=[0.0, 1.0, 0.0]),
        SDFAtom(;coords=[1.41, -0.5, 0.0]),
        SDFAtom(;coords=[-1.41, -0.5, 0.0])
    ]
    bonds = [
        SDFBond(;notation=0),
        SDFBond(;notation=0),
        SDFBond(;notation=0)
    ]
    edges = Edge.([(1,2), (1,3), (1,4)])
    implh_uns = MolGraph(edges, atoms, bonds, on_init=sdf_on_init!)
    @test isempty(implh_uns[:stereocenter])
    @test !haskey(implh_uns[:logs], "warning_stereocenter_ignored")

    bonds = [
        SDFBond(;notation=1),
        SDFBond(;notation=0),
        SDFBond(;notation=0)
    ]
    implh1 = MolGraph(edges, atoms, bonds, on_init=sdf_on_init!)
    @test implh1[:stereocenter][1] == (2, 3, 4, true)

    bonds = [
        SDFBond(;notation=0),
        SDFBond(;notation=6),
        SDFBond(;notation=0)
    ]
    implh2 = MolGraph(edges, atoms, bonds, on_init=sdf_on_init!)
    @test implh2[:stereocenter][1] == (2, 3, 4, false)

    bonds = [
        SDFBond(;notation=6),
        SDFBond(;notation=6),
        SDFBond(;notation=0)
    ]
    implh3 = MolGraph(edges, atoms, bonds, on_init=sdf_on_init!)
    @test implh3[:stereocenter][1] == (2, 3, 4, true)

    bonds = [
        SDFBond(;notation=0),
        SDFBond(;notation=6),
        SDFBond(;notation=1)
    ]
    implh_wrong = MolGraph(edges, atoms, bonds, on_init=sdf_on_init!)
    @test isempty(implh_wrong[:stereocenter])
    @test haskey(implh_wrong[:logs], "warning_stereocenter_ignored")

    bonds = [
        SDFBond(;notation=1),
        SDFBond(;notation=1),
        SDFBond(;notation=1)
    ]
    implh_wrong2 = MolGraph(edges, atoms, bonds, on_init=sdf_on_init!)
    @test isempty(implh_wrong2[:stereocenter])
    @test haskey(implh_wrong2[:logs], "warning_stereocenter_ignored")

    bonds = [
        SDFBond(;notation=6),
        SDFBond(;notation=6),
        SDFBond(;notation=6)
    ]
    implh_wrong3 = MolGraph(edges, atoms, bonds, on_init=sdf_on_init!)
    @test isempty(implh_wrong3[:stereocenter])
    @test haskey(implh_wrong3[:logs], "warning_stereocenter_ignored")

    bonds = [
        SDFBond(;notation=6),
        SDFBond(;notation=1),
        SDFBond(;notation=6)
    ]
    implh_wrong4 = MolGraph(edges, atoms, bonds, on_init=sdf_on_init!)
    @test isempty(implh_wrong4[:stereocenter])
    @test haskey(implh_wrong4[:logs], "warning_stereocenter_ignored")

    # transformed
    atoms = [
        SDFAtom(;coords=[8.9, -6.1, 0.0]),
        SDFAtom(;coords=[8.0, -6.1, 0.0]),
        SDFAtom(;coords=[9.6, -5.7, 0.0]),
        SDFAtom(;coords=[8.5, -6.9, 0.0]),
        SDFAtom(;coords=[9.4, -6.6, 0.0])
    ]
    bonds = [
        SDFBond(;notation=1),
        SDFBond(;notation=0),
        SDFBond(;notation=6),
        SDFBond(;notation=0)
    ]
    edges = Edge.([(1,2), (1,3), (1,4), (1,5)])
    tr1 = MolGraph(edges, atoms, bonds, on_init=sdf_on_init!)
    @test tr1[:stereocenter][1] == (2, 3, 5, true)

    bonds = [
        SDFBond(;notation=6),
        SDFBond(;notation=0),
        SDFBond(;notation=0),
        SDFBond(;notation=6)
    ]
    tr2 = MolGraph(edges, atoms, bonds, on_init=sdf_on_init!)
    @test tr2[:stereocenter][1] == (2, 3, 5, false)
    # global_logger(default_logger)
end

@testset "stereobond_from_sdf2d" begin
    atoms = [
        SDFAtom(;coords=[-0.5, 1.41, 0.0]),
        SDFAtom(;coords=[0.0, 0.0, 0.0]),
        SDFAtom(;coords=[1.0, 0.0, 0.0]),
        SDFAtom(;coords=[1.5, 1.41, 0.0])
    ]
    bonds = [
        SDFBond(;order=1), SDFBond(;order=2), SDFBond(;order=3),
    ]
    edges = Edge.([(1,2), (2,3), (3,4)])
    mol1 = MolGraph(edges, atoms, bonds, on_init=sdf_on_init!)
    @test mol1[:stereobond][Edge(2 => 3)] == (1, 4, true)

    atoms[4] = SDFAtom(;coords=[1.5, -1.41, 0.0])
    mol2 = MolGraph(edges, atoms, bonds, on_init=sdf_on_init!)
    @test mol2[:stereobond][Edge(2 => 3)] == (1, 4, false)

    atoms[4] = SDFAtom(;coords=[2.0, 0.0, 0.0])
    mol3 = MolGraph(edges, atoms, bonds, on_init=sdf_on_init!)
    @test isempty(mol3[:stereobond])
end

@testset "stereobond_from_smiles" begin
    # default_logger = global_logger(ConsoleLogger(stdout, Logging.Debug))
    mol1 = smilestomol("C/C=C\\C")
    @test mol1[:stereobond][Edge(2 => 3)] == (1, 4, true)

    mol2 = smilestomol("C/C=C/C")
    @test mol2[:stereobond][Edge(2 => 3)] == (1, 4, false)

    cis = smilestomol("C(/C)=C/C")
    @test cis[:stereobond][Edge(1 => 3)] == (2, 4, true)  # cis

    conflict = smilestomol("C/C=C(/C)/C")
    @test isempty(conflict[:stereobond])
    @test haskey(conflict[:logs], "warning_stereobond_ignored")
    # serialization check
    @test nv(MolGraph(to_json(conflict))) == 5

    mol3 = smilestomol("C\\C([H])=C([H])/C")
    @test mol3[:stereobond][Edge(2 => 4)] == (1, 6, true)

    mol4 = smilestomol("C/C=C(\\C=C)/C=C(/C)C")
    @test mol4[:stereobond][Edge(2 => 3)] == (1, 4, true)
    @test !haskey(mol4[:stereobond], Edge(4 => 5))
    @test mol4[:stereobond][Edge(6 => 7)] == (3, 8, false)

    # https://github.com/mojaie/MolecularGraph.jl/issues/91
    mol5 = smilestomol("C=C/C(=C/Cl)/Cl")
    @test !haskey(mol5[:stereobond], Edge(1 => 2))
    @test mol5[:stereobond][Edge(3 => 4)] == (2, 5, false)

    # Note: coordgenlib error : CHMMol-> sketcher stereo error: wrong number for RSpriorities
    thiostrepton = smilestomol(raw"C[C@@H](O)[C@@](C)(O)[C@@H](C2=NC6=CS2)NC(C1CSC(/C(NC([C@@H](NC(C3=CS[C@]([C@]4(NC([C@@H](NC(C(NC([C@H](C)NC%10=O)=O)=C)=O)C)=O)C(C5=CSC([C@H]([C@@H](C)OC(C8=CC([C@H](C)O)=C(C=CC(N[C@]([H])%10[C@H](CC)C)[C@@H]9O)C9=N8)=O)NC6=O)=N5)N=C(C7=NC(C(NC(C(NC(C(N)=O)=C)=O)=C)=O)=CS7)CC4)=N3)=O)[C@@H](C)O)=O)=C/C)=N1)=O")
    @test thiostrepton[:stereobond][Edge(21 => 123)] == (20, 124, false)

    linoleic_acid = smilestomol(raw"OC(=O)CCCCCCC/C=C\C/C=C\CCCCC")
    @test linoleic_acid[:stereobond][Edge(11 => 12)] == (10, 13, true)
    @test linoleic_acid[:stereobond][Edge(14 => 15)] == (13, 16, true)
    subg, vmap = induced_subgraph(linoleic_acid, collect(8:20))
    # revmap 8 => 1, 9 => 2, 10 => 3, ...
    @test subg[:stereobond][Edge(4 => 5)] == (3, 6, true)
    @test subg[:stereobond][Edge(7 => 8)] == (6, 9, true)
    rem_vertices!(linoleic_acid, collect(4:12))
    # revmap 1 => 1, 2 => 2, 3 => 3, 20 => 4, 19 => 5, 18 => 6, ...
    @test linoleic_acid[:stereobond][Edge(9 => 10)] == (11, 8, true)

    # global_logger(default_logger)
end

end # stereo
