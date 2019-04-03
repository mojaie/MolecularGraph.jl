
using MolecularGraph:
    forward!, chain!, fragment!, componentquery!


@testset "smarts.molecule" begin

@testset "chain" begin
    state = SmartsParser("CCCCC", false)
    state.node = 1
    push!(state.nodeattrs, SmartsAtom(:Symbol => :C))
    forward!(state)
    chain!(state)
    @test length(state.edges) == 4
    @test length(state.nodeattrs) == 5
    @test length(state.edgeattrs) == 4

    state = SmartsParser("C1CCCCC1", false)
    state.node = 1
    push!(state.nodeattrs, SmartsAtom(:Symbol => :C))
    forward!(state)
    chain!(state)
    @test state.edges[6] == (6, 1)
end

@testset "fragment" begin
    branched1 = SmartsParser("C(C)(C)(C)C", false)
    fragment!(branched1)
    @test branched1.edges[4] == (1, 5)

    branched2 = SmartsParser("CC(C)C(C)C(C)C", false)
    fragment!(branched2)
    @test branched2.edges[7] == (6, 8)

    nested1 = SmartsParser("C(C(C(C)C)C)C", false)
    fragment!(nested1)
    @test nested1.edges[5] == (2, 6)
    @test nested1.edges[6] == (1, 7)

    ring1 = SmartsParser("CC(C(C)1)CC1", false)
    fragment!(ring1)
    @test ring1.edges[6] == (6, 3)

    invalid1 = SmartsParser("CC(", false)
    @test_throws ErrorException fragment!(invalid1)

    invalid2 = SmartsParser("C()C", false)
    @test_throws ErrorException fragment!(invalid2)

    invalid3 = SmartsParser("C1CC", false)
    @test_throws ErrorException fragment!(invalid3)

    invalid4 = SmartsParser("1CCC1", false)
    @test_throws ErrorException fragment!(invalid4)

    invalid5 = SmartsParser("CC(C)", false)
    @test_throws ErrorException fragment!(invalid5)

    invalid6 = SmartsParser("(CC)CC", false)
    @test_throws ErrorException fragment!(invalid6)

    invalid7 = SmartsParser("C(C(C))CC", false)
    @test_throws ErrorException fragment!(invalid7)

    invalid8 = SmartsParser("C(1C)C1C", false)
    @test_throws ErrorException fragment!(invalid8)
end


@testset "component" begin
    conn0 = SmartsParser("C.CC.CCC", true)
    componentquery!(conn0)
    @test isempty(conn0.connectivity)

    conn1 = SmartsParser("(C.CC.CCC)", true)
    componentquery!(conn1)
    @test conn1.connectivity[1] == [1, 2, 4]

    conn2 = SmartsParser("CC.CC.(CC)", true)
    componentquery!(conn2)
    @test conn2.connectivity[1] == [5]

    conn3 = SmartsParser("CC.(CC(C)C.CC)", true)
    componentquery!(conn3)
    @test conn3.connectivity[1] == [3, 7]

    conn4 = SmartsParser("(CCC).(CCC)", true)
    componentquery!(conn4)
    @test conn4.connectivity == [[1], [4]]

    conn5 = SmartsParser("(C.C.C).(C.C).C.(C)", true)
    componentquery!(conn5)
    @test conn5.connectivity == [[1, 2, 3], [4, 5], [7]]

    invalid1 = SmartsParser("C..C", true)
    @test_throws ErrorException fragment!(invalid1)

    invalid2 = SmartsParser("CCC.", true)
    @test_throws ErrorException fragment!(invalid2)

    invalid3 = SmartsParser("CC(C).C", true)
    @test_throws ErrorException fragment!(invalid3)
end

end # smarts.molecule
