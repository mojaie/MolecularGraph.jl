
@testset "smarts.molecule" begin

@testset "chain" begin
    pen = SmartsParserState("CCCCC")
    chain!(pen, 0, true)
    @test atomcount(pen.mol) == 5
    @test bondcount(pen.mol) == 4

    hex = SmartsParserState("C1CCCCC1")
    chain!(hex, 0, true)
    @test 6 in keys(neighbors(hex.mol, 1))
end

@testset "fragment" begin
    neo = SmartsParserState("C(C)(C)(C)C")
    fragment!(neo, 0, true)
    @test degree(neo.mol, 1) == 4

    nested1 = SmartsParserState("C(C(C(C)C)C)C")
    fragment!(nested1, 0, true)
    @test 7 in keys(neighbors(nested1.mol, 1))
    @test 2 in keys(neighbors(nested1.mol, 6))

    invalid1 = SmartsParserState("CC(C)")
    @test_throws MolParseError fragment!(invalid1, 0, true)

    invalid2 = SmartsParserState("(CC)CC")
    @test_throws MolParseError fragment!(invalid2, 0, true)

    invalid3 = SmartsParserState("C(C(C))CC")
    @test_throws MolParseError fragment!(invalid3, 0, true)
end


@testset "component" begin
    conn0 = SmartsParserState("C.CC.CCC")
    componentquery!(conn0)
    @test degree(conn0.mol, 1) == 0
    @test degree(conn0.mol, 3) == 1
    @test isempty(conn0.mol.connectivity)

    conn1 = SmartsParserState("(C.CC.CCC)")
    componentquery!(conn1)
    @test conn1.mol.connectivity[1] == [1, 2, 4]

    conn2 = SmartsParserState("CC.CC.(CC)")
    componentquery!(conn2)
    @test conn2.mol.connectivity[1] == [5]

    conn3 = SmartsParserState("CC.(CC(C)C.CC)")
    componentquery!(conn3)
    @test conn3.mol.connectivity[1] == [3, 7]

    conn4 = SmartsParserState("(CCC).(CCC)")
    componentquery!(conn4)
    @test conn4.mol.connectivity == [[1], [4]]

    conn5 = SmartsParserState("(C.C.C).(C.C).C.(C)")
    componentquery!(conn5)
    @test conn5.mol.connectivity == [[1, 2, 3], [4, 5], [7]]
end

end # smarts.molecule
