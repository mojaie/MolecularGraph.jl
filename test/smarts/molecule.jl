
@testset "smarts.molecule" begin

@testset "chain" begin
    state = ConnectedSmarts("CCCCC")
    state.node = 1
    updateatom!(state.mol, SmartsAtom(:Symbol => :C), 1)
    forward!(state)
    chain!(state)
    @test atomcount(state.mol) == 5
    @test bondcount(state.mol) == 4

    state = ConnectedSmarts("C1CCCCC1")
    state.node = 1
    updateatom!(state.mol, SmartsAtom(:Symbol => :C), 1)
    forward!(state)
    chain!(state)
    @test 6 in keys(neighbors(state.mol, 1))
end

@testset "fragment" begin
    branched1 = ConnectedSmarts("C(C)(C)(C)C")
    fragment!(branched1)
    @test degree(branched1.mol, 1) == 4

    branched2 = ConnectedSmarts("CC(C)C(C)C(C)C")
    fragment!(branched2)
    @test degree(branched2.mol, 6) == 3

    nested1 = ConnectedSmarts("C(C(C(C)C)C)C")
    fragment!(nested1)
    @test 7 in keys(neighbors(nested1.mol, 1))
    @test 2 in keys(neighbors(nested1.mol, 6))

    ring1 = ConnectedSmarts("CC(C(C)1)CC1")
    fragment!(ring1)
    @test 3 in keys(neighbors(ring1.mol, 6))

    invalid1 = ConnectedSmarts("CC(")
    @test_throws MolParseError fragment!(invalid1)

    invalid2 = ConnectedSmarts("C()C")
    @test_throws MolParseError fragment!(invalid2)

    invalid3 = ConnectedSmarts("C1CC")
    @test_throws MolParseError fragment!(invalid3)

    invalid4 = ConnectedSmarts("1CCC1")
    @test_throws MolParseError fragment!(invalid4)

    invalid5 = ConnectedSmarts("CC(C)")
    @test_throws MolParseError fragment!(invalid5)

    invalid6 = ConnectedSmarts("(CC)CC")
    @test_throws MolParseError fragment!(invalid6)

    invalid7 = ConnectedSmarts("C(C(C))CC")
    @test_throws MolParseError fragment!(invalid7)

    invalid8 = ConnectedSmarts("C(1C)C1C")
    @test_throws MolParseError fragment!(invalid8)
end


@testset "component" begin
    conn0 = DisconnectedSmarts("C.CC.CCC")
    componentquery!(conn0)
    @test degree(conn0.mol, 1) == 0
    @test degree(conn0.mol, 3) == 1
    @test isempty(conn0.mol.connectivity)

    conn1 = DisconnectedSmarts("(C.CC.CCC)")
    componentquery!(conn1)
    @test conn1.mol.connectivity[1] == [1, 2, 4]

    conn2 = DisconnectedSmarts("CC.CC.(CC)")
    componentquery!(conn2)
    @test conn2.mol.connectivity[1] == [5]

    conn3 = DisconnectedSmarts("CC.(CC(C)C.CC)")
    componentquery!(conn3)
    @test conn3.mol.connectivity[1] == [3, 7]

    conn4 = DisconnectedSmarts("(CCC).(CCC)")
    componentquery!(conn4)
    @test conn4.mol.connectivity == [[1], [4]]

    conn5 = DisconnectedSmarts("(C.C.C).(C.C).C.(C)")
    componentquery!(conn5)
    @test conn5.mol.connectivity == [[1, 2, 3], [4, 5], [7]]

    invalid1 = DisconnectedSmarts("C..C")
    @test_throws MolParseError fragment!(invalid1)

    invalid2 = DisconnectedSmarts("CCC.")
    @test_throws MolParseError fragment!(invalid2)

    invalid3 = DisconnectedSmarts("CC(C).C")
    @test_throws MolParseError fragment!(invalid3)
end

end # smarts.molecule
