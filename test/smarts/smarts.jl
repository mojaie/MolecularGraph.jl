
@testset "smarts.smarts" begin

@testset "chain" begin
    nullmol = parse(SMARTS, "")
    @test nodecount(nullmol) == 0
    @test edgecount(nullmol) == 0

    aliphatic = parse(SMARTS, "C")
    @test nodeattr(aliphatic, 1).query == (
        :and => (:atomsymbol => :C, :isaromatic => false)
    )
    @test edgecount(aliphatic) == 0

    carbonyl = parse(SMARTS, "[CX3]=[OX1]")
    @test nodeattr(carbonyl, 1).query == (:and => (
        :and => (:atomsymbol => :C, :isaromatic => false),
        :connectivity => 3
    ))
    @test nodeattr(carbonyl, 2).query == (:and => (
        :and => (:atomsymbol => :O, :isaromatic => false),
        :connectivity => 1
    ))

    ether = parse(SMARTS, "[#6][OD2][#6]")
    @test nodeattr(ether, 1).query == (:atomsymbol => :C)
    @test nodeattr(ether, 2).query == (:and => (
        :and => (:atomsymbol => :O, :isaromatic => false),
        :nodedegree => 2
    ))

    notH = parse(SMARTS, "[!#1]")
    @test nodeattr(notH, 1).query == (:not => (:atomsymbol => :H))
end

end # smarts.SMARTS
