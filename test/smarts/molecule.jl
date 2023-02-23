
using MolecularGraph:
    forward!, chain!, fragment!, componentquery!


@testset "smarts.molecule" begin

@testset "chain" begin
    state = SMARTSParser{SMARTSMolGraph}("CCCCC")
    state.node = 1
    push!(state.vprops, QueryTruthTable(v -> v[1], [(:symbol, :C)]))
    forward!(state)
    chain!(state)
    @test length(state.edges) == 4
    @test length(state.vprops) == 5
    @test length(state.eprops) == 4

    state = SMARTSParser{SMARTSMolGraph}("C1CCCCC1")
    state.node = 1
    push!(state.vprops, QueryTruthTable(v -> v[1], [(:symbol, :C)]))
    forward!(state)
    chain!(state)
    @test state.edges[6] == Edge(1, 6)

    state = SMARTSParser{SMARTSMolGraph}("C%10CCCCC%10")
    state.node = 1
    push!(state.vprops, QueryTruthTable(v -> v[1], [(:symbol, :C)]))
    forward!(state)
    chain!(state)
    @test state.edges[6] == Edge(1, 6)
end

@testset "fragment" begin
    branched1 = SMARTSParser{SMARTSMolGraph}("C(C)(C)(C)C")
    fragment!(branched1)
    @test branched1.edges[4] == Edge(1, 5)

    branched2 = SMARTSParser{SMARTSMolGraph}("CC(C)C(C)C(C)C")
    fragment!(branched2)
    @test branched2.edges[7] == Edge(6, 8)

    nested1 = SMARTSParser{SMARTSMolGraph}("C(C(C(C)C)C)C")
    fragment!(nested1)
    @test nested1.edges[5] == Edge(2, 6)
    @test nested1.edges[6] == Edge(1, 7)

    ring1 = SMARTSParser{SMARTSMolGraph}("CC(C(C)1)CC1")
    fragment!(ring1)
    @test ring1.edges[6] == Edge(3, 6)

    invalid1 = SMARTSParser{SMARTSMolGraph}("CC(")
    @test_throws ErrorException fragment!(invalid1)

    invalid2 = SMARTSParser{SMARTSMolGraph}("C()C")
    @test_throws ErrorException fragment!(invalid2)

    invalid3 = SMARTSParser{SMARTSMolGraph}("C1CC")
    @test_throws ErrorException fragment!(invalid3)

    invalid4 = SMARTSParser{SMARTSMolGraph}("1CCC1")
    @test_throws ErrorException fragment!(invalid4)

    valid5 = SMARTSParser{SMARTSMolGraph}("CC(C)")
    fragment!(valid5)

    invalid6 = SMARTSParser{SMARTSMolGraph}("(CC)CC")
    @test_throws ErrorException fragment!(invalid6)

    valid7 = SMARTSParser{SMARTSMolGraph}("C(C(C))CC")
    fragment!(valid7)

    invalid8 = SMARTSParser{SMARTSMolGraph}("C(1C)C1C")
    @test_throws ErrorException fragment!(invalid8)

    valid9 = SMARTSParser{SMARTSMolGraph}("C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C(C))))))))))))))))))))C")
    fragment!(valid9)
    @test length(valid9.vprops) == 22
    @test length(valid9.eprops) == 21

end


@testset "component" begin
    conn0 = SMARTSParser{SMARTSMolGraph}("C.CC.CCC")
    componentquery!(conn0)
    @test isempty(conn0.connectivity)

    conn1 = SMARTSParser{SMARTSMolGraph}("(C.CC.CCC)")
    componentquery!(conn1)
    @test conn1.connectivity[1] == [1, 2, 4]

    conn2 = SMARTSParser{SMARTSMolGraph}("CC.CC.(CC)")
    componentquery!(conn2)
    @test conn2.connectivity[1] == [5]

    conn3 = SMARTSParser{SMARTSMolGraph}("CC.(CC(C)C.CC)")
    componentquery!(conn3)
    @test conn3.connectivity[1] == [3, 7]

    conn4 = SMARTSParser{SMARTSMolGraph}("(CCC).(CCC)")
    componentquery!(conn4)
    @test conn4.connectivity == [[1], [4]]

    conn5 = SMARTSParser{SMARTSMolGraph}("(C.C.C).(C.C).C.(C)")
    componentquery!(conn5)
    @test conn5.connectivity == [[1, 2, 3], [4, 5], [7]]

    invalid1 = SMARTSParser{SMARTSMolGraph}("C..C")
    @test_throws ErrorException fragment!(invalid1)

    invalid2 = SMARTSParser{SMARTSMolGraph}("CCC.")
    @test_throws ErrorException fragment!(invalid2)

    valid3 = SMARTSParser{SMARTSMolGraph}("CC(C).C")
    componentquery!(conn0)
    @test isempty(valid3.connectivity)
end

end # smarts.molecule
