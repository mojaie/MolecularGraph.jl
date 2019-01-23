
@testset "annotation.wclogp" begin

@testset "wclogp" begin
    # Wildman and Crippen 1999, Table.2
    mol1 = smilestomol("C=1C=CC=C(OC)C=1O")
    wclogpcalc!(mol1)
    @test mol1.v[:WCLogP] == [
        :C18, :C18, :C18, :C18, :C23, :O4, :C3, :C23, :O2
    ]
    @test wclogp(mol1) == 1.40

    mol2 = smilestomol("C1=CC=CC=C1C2=CC=CC=N2")
    wclogpcalc!(mol2)
    @test mol2.v[:WCLogP] == [
        :C18, :C18, :C18, :C18, :C18, :C20, :C20,
        :C18, :C18, :C18, :C18, :N11
    ]
    @test wclogp(mol2) == 2.75

    # Test molecules
    TESTMOL_DIR = joinpath(dirname(@__FILE__), "..", "..", "assets", "test")
    aliphatic = sdftomol(open(joinpath(TESTMOL_DIR, "wctype_C_alip.mol")))
    wclogpcalc!(aliphatic)
    @test aliphatic.v[:WCLogP] == [
        :C2, :C1, :C7, :C4, :N2, :C3, :N1, :C7, :C5, :O9,
        :C27, :Me1, :C27, :C27, :C27, :H1, :C6, :C26, :C21, :C21,
        :N12, :C18, :C21, :C21, :C11, :O3, :C12, :C1, :C1, :C1,
        :C9, :C10, :C1
    ]

    aromatic = sdftomol(open(joinpath(TESTMOL_DIR, "wctype_C_arom.mol")))
    wclogpcalc!(aromatic)
    @test aromatic.v[:WCLogP] == [
        :C21, :C23, :C22, :C24, :C19, :C19, :C14, :C20, :C25, :O1,
        :O8, :O2, :N3, :S1, :C20, :C15, :C16, :C17, :C13, :C18,
        :Cl, :Br, :I, :C8, :P, :F
    ]

    nitrogen = sdftomol(open(joinpath(TESTMOL_DIR, "wctype_N.mol")))
    wclogpcalc!(nitrogen)
    @test nitrogen.v[:WCLogP] == [
        :N11, :C21, :C22, :C18, :C22, :C22, :N4, :C4, :N8, :C4,
        :N6, :N7, :N5, :C3, :N13, :C3, :C3, :C3, :C3, :N6,
        :N10, :C7, :N9, :N14, :N14
    ]

    oands = sdftomol(open(joinpath(TESTMOL_DIR, "wctype_OS.mol")))
    wclogpcalc!(oands)
    @test oands.v[:WCLogP] == [
        :S2, :C4, :N13, :P, :C4, :O6, :O5, :O7, :C5, :O9,
        :O3, :C5, :O4, :C23, :S3, :C23, :C24, :C21, :O2, :O11,
        :C5, :O12, :O10, :O4, :S1, :O2, :O2
    ]
end

end # annotation.wclogp
