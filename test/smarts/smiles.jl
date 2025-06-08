
@testset "smarts.smiles" begin

@testset "chain" begin
    nullmol = smilestomol("")
    @test nv(nullmol) == 0
    @test ne(nullmol) == 0

    methane = smilestomol("C")
    @test get_prop(methane, 1, :symbol) === :C
    @test ne(methane) == 0

    EtOH = smilestomol("CCO")
    @test get_prop(EtOH, 3, :symbol) === :O

    EtOH2 = smilestomol("[CH3][CH2][OH]")
    @test nv(EtOH2) == 9
    @test ne(EtOH2) == 8

    nitro = smilestomol("N#N")
    @test get_prop(nitro, 1, 2, :order) == 3

    vinylcl = smilestomol("C=CCl")
    @test get_prop(vinylcl, 1, 2, :order) == 2
    @test get_prop(vinylcl, 3, :symbol) === :Cl
end

@testset "branch" begin
    tBuOH = smilestomol("CC(C)(C)O")
    @test get_prop(tBuOH, 5, :symbol) === :O
    @test degree(tBuOH, 2) == 4

    TFA = smilestomol("FC(F)(F)C(=O)O")
    @test get_prop(TFA, 5, 6, :order) == 2
    @test degree(TFA, 2) == 4

    EDTA = smilestomol("C(CN(CC(=O)O)CC(=O)O)N(CC(=O)O)CC(=O)O")
    @test degree(EDTA, 5) == 3
    @test degree(EDTA, 12) == 3
    @test degree(EDTA, 14) == 3
end

@testset "ring" begin
    cychexe = smilestomol("C1C=CCCC1")
    @test degree(cychexe, 1) == 2

    cychexe2 = smilestomol("C=1CCCCC1")
    @test get_prop(cychexe2, 6, 1, :order) == 2

    cychexe3 = smilestomol("C1CCCCC=1")
    @test get_prop(cychexe3, 6, 1, :order) == 2

    cychexe4 = smilestomol("C=1CCCCC=1")
    @test get_prop(cychexe4, 6, 1, :order) == 2

    bicyclohexyl = smilestomol("C1CCCCC1C1CCCCC1")
    @test get_prop(bicyclohexyl, 6, 1, :order) == 1
    @test get_prop(bicyclohexyl, 12, 7, :order) == 1

    bicyclohexyl2 = smilestomol("C1CCCCC1C2CCCCC2")
    @test get_prop(bicyclohexyl2, 6, 1, :order) == 1
    @test get_prop(bicyclohexyl2, 12, 7, :order) == 1

    naphthalene = smilestomol("C=1C=CC=C2C1C=CC=C2")
    @test get_prop(naphthalene, 6, 1, :order) == 2
    @test get_prop(naphthalene, 10, 5, :order) == 1

    naphthalene2 = smilestomol("C1=CC=C2C(=C1)C=CC=C2")
    @test get_prop(naphthalene2, 6, 1, :order) == 1
    @test get_prop(naphthalene2, 10, 4, :order) == 1

    tetrahedrane = smilestomol("C12C3C1C23")
    @test ne(tetrahedrane) == 6

    thiamin = smilestomol("OCCc1c(C)[n+](=cs1)Cc2cnc(C)nc(N)2")
    @test ne(thiamin) == 19

    c60 = smilestomol("c12c3c4c5c1c6c7c8c2c9c1c3c2c3c4c4c%10c5c5c6c6c7c7c%11c8c9c8c9c1c2c1c2c3c4c3c4c%10c5c5c6c6c7c7c%11c8c8c9c1c1c2c3c2c4c5c6c3c7c8c1c23")
    @test ne(c60) == 90
end

@testset "component" begin
    NaCl = smilestomol("[Na+].[Cl-]")
    @test get_prop(NaCl, 1, :symbol) === :Na
    @test get_prop(NaCl, 1, :charge) == 1
    @test get_prop(NaCl, 2, :symbol) === :Cl
    @test get_prop(NaCl, 2, :charge) == -1

    ferrous = smilestomol("[Fe++]")
    @test get_prop(ferrous, 1, :symbol) === :Fe
    @test get_prop(ferrous, 1, :charge) == 2

    ferrous2 = smilestomol("[Fe+2]")
    @test get_prop(ferrous2, 1, :symbol) === :Fe
    @test get_prop(ferrous2, 1, :charge) == 2

    CuSO4 = smilestomol("[Cu+2].[O-]S(=O)(=O)[O-]")
    @test get_prop(CuSO4, 1, :charge) == 2
    @test get_prop(CuSO4, 6, :charge) == -1
    @test degree(CuSO4, 3) == 4

    PhONa = smilestomol("c1cc([O-].[Na+])ccc1")
    @test degree(PhONa, 3) == 3
    @test degree(PhONa, 5) == 0
end

@testset "stereo" begin
    LAla = smilestomol("N[C@H](C)C(=O)O")
    @test degree(LAla, 2) == 4
    @test get_prop(LAla, 2, :stereo) === :anticlockwise
    DAla = smilestomol("N[C@@H](C)C(=O)O")
    @test degree(DAla, 2) == 4
    @test get_prop(DAla, 2, :stereo) === :clockwise
end

@testset "aromatic" begin
    benzene = smilestomol("c1ccccc1")
    @test get_prop(benzene, 1, :isaromatic)
    @test get_prop(benzene, 1, 6, :order) == 1
    pyrrole = smilestomol("[nH]1cccc1")
    @test get_prop(pyrrole, 1, :isaromatic)
    @test degree(pyrrole, 1) == 3
end

end # smarts.smiles
