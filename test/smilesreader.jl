
@testset "smilesreader" begin

@testset "tokenize" begin
    @test tokenize("") == []
    @test tokenize("CC(=O)CC#N") == ["C", "C", "(", "=O", ")", "C", "C", "#N"]
    @test tokenize("[CH3][CH2][O-].[Na+]") == [
        "[CH3]", "[CH2]", "[O-]", ".", "[Na+]"]
    @test tokenize("C=1CCCCC=1") == ["C=1", "C", "C", "C", "C", "C=1"]
    @test tokenize("c1nopc1") == ["c1", "n", "o", "p", "c1"]
    @test tokenize("C\\C=C/C") == ["C", "\\C", "=C", "/C"]
    @test tokenize("BrCCl") == ["Br", "C", "Cl"]
    @test tokenize("C12C3C1C23") == ["C12", "C3", "C1", "C23"]
    @test_throws OperationError tokenize("C1CCC(C)1CC")
    @test_throws OperationError tokenize("MgCl")
end

@testset "parsetoken" begin
    @test parsesmilestoken("=O") == ("O", "=", nothing)
    @test parsesmilestoken("[13C]=1") == ("13C", nothing, "=1")
    @test parsesmilestoken("#n2") == ("n", "#", "2")
    @test parsesmilestoken("P=1=2") == ("P", nothing, "=1=2")
end

@testset "parseatom" begin
    br = parsesmilesatom("Br")
    @test br[1] == :Br
    iso = parsesmilesatom("14cH2")
    @test iso[1] == :C
    @test iso[4] == 14.0
    @test iso[5]
    zn = parsesmilesatom("Zn++")
    @test zn[1] == :Zn
    @test zn[2] == 2
    ani = parsesmilesatom("O-2")
    @test ani[2] == -2
end

@testset "parsebond" begin
    tri = parsesmilesbond("#")
    @test tri[1] == 3
end

@testset "smilestomol" begin
    nullmol = smilestomol("")
    @test atomcount(nullmol) == 0
    @test bondcount(nullmol) == 0
    methane = smilestomol("C")
    @test getatom(methane, 1).symbol == :C
    @test bondcount(methane) == 0
    EtOH = smilestomol("CCO")
    @test getatom(EtOH, 3).symbol == :O
    EtOH2 = smilestomol("[CH3][CH2][OH]")
    @test getbond(EtOH2, 2, 3).order == 1
    nitro = smilestomol("N#N")
    @test getbond(nitro, 1, 2).order == 3
    vinylcl = smilestomol("C=CCl")
    @test getbond(vinylcl, 1, 2).order == 2
    @test getatom(vinylcl, 3).symbol == :Cl
    tetrahedrane = smilestomol("C12C3C1C23")
    @test bondcount(tetrahedrane) == 6

    tBuOH = smilestomol("CC(C)(C)O")
    @test getatom(tBuOH, 5).symbol == :O
    @test neighborcount(tBuOH, 2) == 4
    TFA = smilestomol("FC(F)(F)C(=O)O")
    @test getbond(TFA, 5, 6).order == 2
    @test neighborcount(TFA, 2) == 4
    EDTA = smilestomol("C(CN(CC(=O)O)CC(=O)O)N(CC(=O)O)CC(=O)O")
    @test neighborcount(EDTA, 5) == 3
    @test neighborcount(EDTA, 12) == 3
    @test neighborcount(EDTA, 14) == 3

    cychexe = smilestomol("C1C=CCCC1")
    @test neighborcount(cychexe, 1) == 2
    cychexe2 = smilestomol("C=1CCCCC1")
    @test getbond(cychexe2, 1, 6).order == 2
    cychexe3 = smilestomol("C1CCCCC=1")
    @test getbond(cychexe3, 1, 6).order == 2
    bicyclohexyl = smilestomol("C1CCCCC1C1CCCCC1")
    @test getbond(bicyclohexyl, 1, 6).order == 1
    @test getbond(bicyclohexyl, 7, 12).order == 1
    bicyclohexyl2 = smilestomol("C1CCCCC1C2CCCCC2")
    @test getbond(bicyclohexyl2, 1, 6).order == 1
    @test getbond(bicyclohexyl2, 7, 12).order == 1
    naphthalene = smilestomol("C=1C=CC=C2C1C=CC=C2")
    @test getbond(naphthalene, 1, 6).order == 2
    @test getbond(naphthalene, 5, 10).order == 1
    naphthalene2 = smilestomol("C1=CC=C2C(=C1)C=CC=C2")
    @test getbond(naphthalene2, 1, 6).order == 1
    @test getbond(naphthalene2, 4, 10).order == 1

    NaCl = smilestomol("[Na+].[Cl-]")
    @test getatom(NaCl, 1).symbol == :Na
    @test getatom(NaCl, 1).charge == 1
    @test getatom(NaCl, 2).symbol == :Cl
    @test getatom(NaCl, 2).charge == -1
    ferrous = smilestomol("[Fe++]")
    @test getatom(ferrous, 1).symbol == :Fe
    @test getatom(ferrous, 1).charge == 2
    ferrous2 = smilestomol("[Fe+2]")
    @test getatom(ferrous2, 1).symbol == :Fe
    @test getatom(ferrous2, 1).charge == 2
    CuSO4 = smilestomol("[Cu+2].[O-]S(=O)(=O)[O-]")
    @test getatom(CuSO4, 1).charge == 2
    @test getatom(CuSO4, 6).charge == -1
    @test neighborcount(CuSO4, 3) == 4

    # TODO: optical isomers
    LAla = smilestomol("N[C@H](C)C(=O)O")
    @test neighborcount(LAla, 2) == 3
    @test getatom(LAla, 2).smiles_stereo == 1
    DAla = smilestomol("N[C@@H](C)C(=O)O")
    @test neighborcount(DAla, 2) == 3
    @test getatom(DAla, 2).smiles_stereo == 2

    # TODO: aromatic
    benzene = smilestomol("c1ccccc1")
    @test getatom(benzene, 1).smiles_aromatic
    @test getbond(benzene, 1, 6).order == 1
    pyridone = smilestomol("O=c1[nH]cccc1")
    @test getatom(pyridone, 3).smiles_aromatic
    @test neighborcount(pyridone, 3) == 2 # Implicit H is ignored
    PhONa = smilestomol("c1cc([O-].[Na+])ccc1")
    @test neighborcount(PhONa, 3) == 3
    @test neighborcount(PhONa, 5) == 0

    # println(Int64[a.index for a in atomvector(ethanol2)])
    # println(Tuple{Int64, Int64}[(b.u, b.v) for b in bondvector(ethanol2)])
end

end # smilesreader
