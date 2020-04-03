
@testset "stereo" begin

@testset "stereohydrogens" begin
    LAla1 = parse(SMILES, "N[C@@H](C)C(=O)O")
    @test LAla1.nodeattrs[2].stereo == :clockwise
    LAla1 = removestereohydrogens(LAla1)
    @test nodecount(LAla1) == 6
    @test LAla1.nodeattrs[2].stereo == :clockwise
    @test LAla1.nodeattrs[3].symbol == :C
    LAla1 = addstereohydrogens(LAla1)
    @test nodecount(LAla1) == 7
    @test LAla1.nodeattrs[2].stereo == :clockwise
    @test LAla1.nodeattrs[3].symbol == :H

    LAla2 = parse(SMILES, "N[C@@](C)([H])C(=O)O")
    LAla2 = removestereohydrogens(LAla2)
    @test nodecount(LAla2) == 6
    @test LAla2.nodeattrs[2].stereo == :anticlockwise
    
    LAla3 = parse(SMILES, "[H][C@@](C)(N)C(=O)O")
    LAla3 = removestereohydrogens(LAla3)
    @test nodecount(LAla3) == 6
    @test LAla3.nodeattrs[1].stereo == :clockwise

    LAla4 = parse(SMILES, "N[C@@](C)(C(=O)O)[H]")
    LAla4 = removestereohydrogens(LAla4)
    @test nodecount(LAla4) == 6
    @test LAla4.nodeattrs[2].stereo == :clockwise
end

@testset "setdiastereosmiles" begin
    mol1 = parse(SMILES, "C/C=C\\C")
    setdiastereo!(mol1)
    @test mol1.edgeattrs[2].stereo === :cis

    mol2 = parse(SMILES, "C/C=C/C")
    setdiastereo!(mol2)
    @test mol2.edgeattrs[2].stereo === :trans

    mol3 = parse(SMILES, "C\\C([H])=C([H])/C")
    setdiastereo!(mol3)
    @test mol3.edgeattrs[3].stereo === :trans

    mol4 = parse(SMILES, "C/C=C(\\C=C)/C=C(/C)C")
    setdiastereo!(mol4)
    @test mol4.edgeattrs[2].stereo === :cis
    @test mol4.edgeattrs[4].stereo === :unspecified
    @test mol4.edgeattrs[6].stereo === :trans
end

@testset "setdiastereosdfile" begin
    bonds = [
        SDFileBond(1), SDFileBond(2), SDFileBond(3),
    ]
    
    atoms = [
        SDFileAtom(:C, 0, 1, nothing, [-0.5, 1.41]),
        SDFileAtom(:C, 0, 1, nothing, [0.0, 0.0]),
        SDFileAtom(:C, 0, 1, nothing, [1.0, 0.0]),
        SDFileAtom(:C, 0, 1, nothing, [1.5, 1.41])
    ]
    mol1 = graphmol([(1,2), (2,3), (3,4)], atoms, bonds)
    setdiastereo!(mol1)
    @test mol1.edgeattrs[2].stereo === :cis

    atoms[4] = SDFileAtom(:C, 0, 1, nothing, [1.5, -1.41])
    mol2 = graphmol([(1,2), (2,3), (3,4)], atoms, bonds)
    setdiastereo!(mol2)
    @test mol2.edgeattrs[2].stereo === :trans

    atoms[4] = SDFileAtom(:C, 0, 1, nothing, [2.0, 0.0])
    mol3 = graphmol([(1,2), (2,3), (3,4)], atoms, bonds)
    setdiastereo!(mol3)
    @test mol3.edgeattrs[2].stereo === :unspecified
end

@testset "setstereocenter" begin
    atoms = [
        SDFileAtom(:C, 0, 1, nothing, [0.0, 0.0]),
        SDFileAtom(:C, 0, 1, nothing, [-1.41, 0.5]),
        SDFileAtom(:C, 0, 1, nothing, [1.41, 0.5]),
        SDFileAtom(:C, 0, 1, nothing, [-0.5, -1.41]),
        SDFileAtom(:C, 0, 1, nothing, [0.5, -1.41])
    ]
    bonds = [
        SDFileBond(1, 0), SDFileBond(1, 0), SDFileBond(1, 1), SDFileBond(1, 6)
    ]
    mol1 = graphmol([(1,2), (1,3), (1,4), (1,5)], atoms, bonds)
    setstereocenter!(mol1)
    @test mol1.nodeattrs[1].stereo === :clockwise

    bonds = [
        SDFileBond(1, 1), SDFileBond(1, 0), SDFileBond(1, 6), SDFileBond(1, 0)
    ]
    mol2 = graphmol([(1,2), (1,3), (1,4), (1,5)], atoms, bonds)
    setstereocenter!(mol2)
    @test mol2.nodeattrs[1].stereo === :anticlockwise

    bonds = [
        SDFileBond(1, 1), SDFileBond(1, 0), SDFileBond(1, 0), SDFileBond(1, 1)
    ]
    mol3 = graphmol([(1,2), (1,3), (1,4), (1,5)], atoms, bonds)
    setstereocenter!(mol3)
    @test mol3.nodeattrs[1].stereo === :anticlockwise

    bonds = [
        SDFileBond(1, 6), SDFileBond(1, 0), SDFileBond(1, 0), SDFileBond(1, 6)
    ]
    mol4 = graphmol([(1,2), (1,3), (1,4), (1,5)], atoms, bonds)
    setstereocenter!(mol4)
    @test mol4.nodeattrs[1].stereo === :clockwise

    bonds = [
        SDFileBond(1, 0), SDFileBond(1, 1), SDFileBond(1, 0), SDFileBond(1, 0)
    ]
    mol5 = graphmol([(1,2), (1,3), (1,4), (1,5)], atoms, bonds)
    setstereocenter!(mol5)
    @test mol5.nodeattrs[1].stereo === :clockwise

    bonds = [
        SDFileBond(1, 0), SDFileBond(1, 0), SDFileBond(1, 6), SDFileBond(1, 0)
    ]
    mol6 = graphmol([(1,2), (1,3), (1,4), (1,5)], atoms, bonds)
    setstereocenter!(mol6)
    @test mol6.nodeattrs[1].stereo === :anticlockwise
end

end # stereo
