
@testset "smarts.molecule" begin

@testset "chain" begin
    state = SmartsParser{SMARTS}("CCCCC", false)
    state.node = 1
    updateatom!(state.mol, SmartsAtom(:Symbol => :C), 1)
    forward!(state)
    chain!(state)
    @test atomcount(state.mol) == 5
    @test bondcount(state.mol) == 4

    state = SmartsParser{SMARTS}("C1CCCCC1", false)
    state.node = 1
    updateatom!(state.mol, SmartsAtom(:Symbol => :C), 1)
    forward!(state)
    chain!(state)
    @test 6 in keys(neighbors(state.mol, 1))
end

@testset "fragment" begin
    branched1 = SmartsParser{SMARTS}("C(C)(C)(C)C", false)
    fragment!(branched1)
    @test degree(branched1.mol, 1) == 4

    branched2 = SmartsParser{SMARTS}("CC(C)C(C)C(C)C", false)
    fragment!(branched2)
    @test degree(branched2.mol, 6) == 3

    nested1 = SmartsParser{SMARTS}("C(C(C(C)C)C)C", false)
    fragment!(nested1)
    @test 7 in keys(neighbors(nested1.mol, 1))
    @test 2 in keys(neighbors(nested1.mol, 6))

    ring1 = SmartsParser{SMARTS}("CC(C(C)1)CC1", false)
    fragment!(ring1)
    @test 3 in keys(neighbors(ring1.mol, 6))

    invalid1 = SmartsParser{SMARTS}("CC(", false)
    @test_throws ErrorException fragment!(invalid1)

    invalid2 = SmartsParser{SMARTS}("C()C", false)
    @test_throws ErrorException fragment!(invalid2)

    invalid3 = SmartsParser{SMARTS}("C1CC", false)
    @test_throws ErrorException fragment!(invalid3)

    invalid4 = SmartsParser{SMARTS}("1CCC1", false)
    @test_throws ErrorException fragment!(invalid4)

    invalid5 = SmartsParser{SMARTS}("CC(C)", false)
    @test_throws ErrorException fragment!(invalid5)

    invalid6 = SmartsParser{SMARTS}("(CC)CC", false)
    @test_throws ErrorException fragment!(invalid6)

    invalid7 = SmartsParser{SMARTS}("C(C(C))CC", false)
    @test_throws ErrorException fragment!(invalid7)

    invalid8 = SmartsParser{SMARTS}("C(1C)C1C", false)
    @test_throws ErrorException fragment!(invalid8)
end


@testset "component" begin
    conn0 = SmartsParser{SMARTS}("C.CC.CCC", true)
    componentquery!(conn0)
    @test degree(conn0.mol, 1) == 0
    @test degree(conn0.mol, 3) == 1
    @test isempty(conn0.mol.connectivity)

    conn1 = SmartsParser{SMARTS}("(C.CC.CCC)", true)
    componentquery!(conn1)
    @test conn1.mol.connectivity[1] == [1, 2, 4]

    conn2 = SmartsParser{SMARTS}("CC.CC.(CC)", true)
    componentquery!(conn2)
    @test conn2.mol.connectivity[1] == [5]

    conn3 = SmartsParser{SMARTS}("CC.(CC(C)C.CC)", true)
    componentquery!(conn3)
    @test conn3.mol.connectivity[1] == [3, 7]

    conn4 = SmartsParser{SMARTS}("(CCC).(CCC)", true)
    componentquery!(conn4)
    @test conn4.mol.connectivity == [[1], [4]]

    conn5 = SmartsParser{SMARTS}("(C.C.C).(C.C).C.(C)", true)
    componentquery!(conn5)
    @test conn5.mol.connectivity == [[1, 2, 3], [4, 5], [7]]

    invalid1 = SmartsParser{SMARTS}("C..C", true)
    @test_throws ErrorException fragment!(invalid1)

    invalid2 = SmartsParser{SMARTS}("CCC.", true)
    @test_throws ErrorException fragment!(invalid2)

    invalid3 = SmartsParser{SMARTS}("CC(C).C", true)
    @test_throws ErrorException fragment!(invalid3)
end

end # smarts.molecule
