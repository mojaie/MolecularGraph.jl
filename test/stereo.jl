#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "stereo" begin

@testset "stereo_hydrogen" begin
    LAla1 = smilestomol("N[C@@H](C)C(=O)O")
    @test get_prop(LAla1, 2, :stereo) === :clockwise
    stereocenter_from_smiles!(LAla1)
    @test get_prop(LAla1, :stereocenter)[2] == (1, 3, 4, true)
    @test remove_stereo_hydrogen!(LAla1, 2)
    @test get_prop(LAla1, :stereocenter)[2] == (1, 4, 5, true)

    LAla2 = smilestomol("N[C@](C)([H])C(=O)O")
    @test get_prop(LAla2, 2, :stereo) === :anticlockwise
    stereocenter_from_smiles!(LAla2)
    @test get_prop(LAla2, :stereocenter)[2] == (1, 3, 4, false)
    @test remove_stereo_hydrogen!(LAla2, 2)
    @test get_prop(LAla2, :stereocenter)[2] == (1, 3, 5, true)
    
    LAla3 = smilestomol("[H][C@@](C(=O)O)(C)N")
    @test get_prop(LAla3, 2, :stereo) === :clockwise
    stereocenter_from_smiles!(LAla3)
    @test get_prop(LAla3, :stereocenter)[2] == (1, 3, 6, true)
    @test remove_stereo_hydrogen!(LAla3, 2)
    @test get_prop(LAla3, :stereocenter)[2] == (3, 6, 1, false)

    LAla4 = smilestomol("OC(=O)[C@]([H])(C)N")
    @test get_prop(LAla4, 4, :stereo) === :anticlockwise
    stereocenter_from_smiles!(LAla4)
    @test get_prop(LAla4, :stereocenter)[4] == (2, 5, 6, false)
    @test remove_stereo_hydrogen!(LAla4, 4)
    @test get_prop(LAla4, :stereocenter)[4] == (2, 6, 5, false)
end

@testset "stereocenter_from_sdf2d" begin
    # degree=4
    # default_logger = global_logger(ConsoleLogger(stdout, Logging.Debug))
    atoms = [
        SDFAtom(:C, 0, 1, nothing, [0.0, 0.0]),
        SDFAtom(:C, 0, 1, nothing, [-1.41, 0.5]),
        SDFAtom(:C, 0, 1, nothing, [1.41, 0.5]),
        SDFAtom(:C, 0, 1, nothing, [-0.5, -1.41]),
        SDFAtom(:C, 0, 1, nothing, [0.5, -1.41])
    ]
    bonds = [
        SDFBond(1, 0),
        SDFBond(1, 0),
        SDFBond(1, 0),
        SDFBond(1, 0)
    ]
    edges = Edge.([(1,2), (1,3), (1,4), (1,5)])
    uns1 = MolGraph(edges, atoms, bonds)
    stereocenter_from_sdf2d!(uns1)
    @test !has_prop(uns1, :stereocenter)

    bonds = [
        SDFBond(1, 0),
        SDFBond(1, 0),
        SDFBond(1, 1),
        SDFBond(1, 6)
    ]
    mol1 = MolGraph(edges, atoms, bonds)
    stereocenter_from_sdf2d!(mol1)
    @test get_prop(mol1, :stereocenter)[1] === (2, 3, 5, false)

    bonds = [
        SDFBond(1, 1),
        SDFBond(1, 0),
        SDFBond(1, 6),
        SDFBond(1, 0)
    ]
    mol2 = MolGraph(edges, atoms, bonds)
    stereocenter_from_sdf2d!(mol2)
    @test get_prop(mol2, :stereocenter)[1] === (2, 3, 5, true)

    bonds = [
        SDFBond(1, 1),
        SDFBond(1, 0),
        SDFBond(1, 0),
        SDFBond(1, 1)
    ]
    mol3 = MolGraph(edges, atoms, bonds)
    stereocenter_from_sdf2d!(mol3)
    @test get_prop(mol3, :stereocenter)[1] === (2, 3, 5, true)

    bonds = [
        SDFBond(1, 6),
        SDFBond(1, 0),
        SDFBond(1, 0),
        SDFBond(1, 6)
    ]
    mol4 = MolGraph(edges, atoms, bonds)
    stereocenter_from_sdf2d!(mol4)
    @test get_prop(mol4, :stereocenter)[1] === (2, 3, 5, false)

    bonds = [
        SDFBond(1, 0),
        SDFBond(1, 1),
        SDFBond(1, 0),
        SDFBond(1, 0)
    ]
    mol5 = MolGraph(edges, atoms, bonds)
    stereocenter_from_sdf2d!(mol5)
    @test get_prop(mol5, :stereocenter)[1] === (2, 3, 5, false)

    bonds = [
        SDFBond(1, 0),
        SDFBond(1, 0),
        SDFBond(1, 6),
        SDFBond(1, 0)
    ]
    mol6 = MolGraph(edges, atoms, bonds)
    stereocenter_from_sdf2d!(mol6)
    @test get_prop(mol6, :stereocenter)[1] === (2, 3, 5, true)

    bonds = [
        SDFBond(1, 6),
        SDFBond(1, 6),
        SDFBond(1, 6),
        SDFBond(1, 0)
    ]
    wrong1 = MolGraph(edges, atoms, bonds)
    stereocenter_from_sdf2d!(wrong1)
    @test !has_prop(wrong1, :stereocenter)

    bonds = [
        SDFBond(1, 1),
        SDFBond(1, 0),
        SDFBond(1, 1),
        SDFBond(1, 1)
    ]
    wrong2 = MolGraph(edges, atoms, bonds)
    stereocenter_from_sdf2d!(wrong2)
    @test !has_prop(wrong2, :stereocenter)

    bonds = [
        SDFBond(1, 1),
        SDFBond(1, 0),
        SDFBond(1, 0),
        SDFBond(1, 6)
    ]
    wrong3 = MolGraph(edges, atoms, bonds)
    stereocenter_from_sdf2d!(wrong3)
    @test !has_prop(wrong3, :stereocenter)

    bonds = [
        SDFBond(1, 1),
        SDFBond(1, 1),
        SDFBond(1, 0),
        SDFBond(1, 0)
    ]
    wrong4 = MolGraph(edges, atoms, bonds)
    stereocenter_from_sdf2d!(wrong4)
    @test !has_prop(wrong4, :stereocenter)

    # degree=3
    atoms = [
        SDFAtom(:C, 0, 1, nothing, [0.0, 0.0]),
        SDFAtom(:C, 0, 1, nothing, [0.0, 1.0]),
        SDFAtom(:C, 0, 1, nothing, [1.41, -0.5]),
        SDFAtom(:C, 0, 1, nothing, [-1.41, -0.5])
    ]
    bonds = [
        SDFBond(1, 0),
        SDFBond(1, 0),
        SDFBond(1, 0)
    ]
    edges = Edge.([(1,2), (1,3), (1,4)])
    uns2 = MolGraph(edges, atoms, bonds)
    stereocenter_from_sdf2d!(uns2)
    @test !has_prop(uns2, :stereocenter)

    bonds = [
        SDFBond(1, 1), SDFBond(1, 0), SDFBond(1, 0)
    ]
    implh1 = MolGraph(edges, atoms, bonds)
    stereocenter_from_sdf2d!(implh1)
    @test get_prop(implh1, :stereocenter)[1] === (2, 3, 4, true)

    bonds = [
        SDFBond(1, 0), SDFBond(1, 6), SDFBond(1, 0)
    ]
    implh2 = MolGraph(edges, atoms, bonds)
    stereocenter_from_sdf2d!(implh2)
    @test get_prop(implh2, :stereocenter)[1] === (2, 3, 4, false)

    bonds = [
        SDFBond(1, 0), SDFBond(1, 6), SDFBond(1, 1)
    ]
    wrong5 = MolGraph(edges, atoms, bonds)
    stereocenter_from_sdf2d!(wrong5)
    @test !has_prop(wrong5, :stereocenter)
    # global_logger(default_logger)
end

@testset "stereobond_from_sdf2d" begin
    atoms = [
        SDFAtom(:C, 0, 1, nothing, [-0.5, 1.41]),
        SDFAtom(:C, 0, 1, nothing, [0.0, 0.0]),
        SDFAtom(:C, 0, 1, nothing, [1.0, 0.0]),
        SDFAtom(:C, 0, 1, nothing, [1.5, 1.41])
    ]
    bonds = [
        SDFBond(1), SDFBond(2), SDFBond(3),
    ]
    edges = Edge.([(1,2), (2,3), (3,4)])
    mol1 = MolGraph(edges, atoms, bonds)
    stereobond_from_sdf2d!(mol1)
    @test get_prop(mol1, :stereobond)[Edge(2 => 3)] === (1, 4, true)

    atoms[4] = SDFAtom(:C, 0, 1, nothing, [1.5, -1.41])
    mol2 = MolGraph(edges, atoms, bonds)
    stereobond_from_sdf2d!(mol2)
    @test get_prop(mol2, :stereobond)[Edge(2 => 3)] === (1, 4, false)

    atoms[4] = SDFAtom(:C, 0, 1, nothing, [2.0, 0.0])
    mol3 = MolGraph(edges, atoms, bonds)
    stereobond_from_sdf2d!(mol3)
    @test !has_prop(mol3, :stereobond)
end

@testset "stereobond_from_smiles" begin
    # default_logger = global_logger(ConsoleLogger(stdout, Logging.Debug))
    mol1 = smilestomol("C/C=C\\C")
    stereobond_from_smiles!(mol1)
    @test get_prop(mol1, :stereobond)[Edge(2 => 3)] === (1, 4, true)

    mol2 = smilestomol("C/C=C/C")
    stereobond_from_smiles!(mol2)
    @test get_prop(mol2, :stereobond)[Edge(2 => 3)] === (1, 4, false)

    # follows OpenSMILES specification http://opensmiles.org/opensmiles.html#chirality
    # -> "up-ness" or "down-ness" of each single bond is relative to the carbon atom
    cis = smilestomol("C(/C)=C/C")
    stereobond_from_smiles!(cis)
    @test get_prop(cis, :stereobond)[Edge(1 => 3)] === (2, 4, true)  # cis

    mol3 = smilestomol("C\\C([H])=C([H])/C")
    stereobond_from_smiles!(mol3)
    @test get_prop(mol3, :stereobond)[Edge(2 => 4)] === (1, 6, true)

    mol4 = smilestomol("C/C=C(\\C=C)/C=C(/C)C")
    stereobond_from_smiles!(mol4)
    @test get_prop(mol4, :stereobond)[Edge(2 => 3)] === (1, 4, true)
    @test !haskey(get_prop(mol4, :stereobond), Edge(4 => 5))
    @test get_prop(mol4, :stereobond)[Edge(6 => 7)] === (3, 8, false)
    # global_logger(default_logger)
end

end # stereo
