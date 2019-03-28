
@testset "smarts.smarts" begin

@testset "chain" begin
    nullmol = parse(SMARTS, "")
    @test atomcount(nullmol) == 0
    @test bondcount(nullmol) == 0

    aliphatic = parse(SMARTS, "C")
    @test getatom(aliphatic, 1).query == (
        :and => (:atomsymbol => :C, :isaromatic => false)
    )
    @test bondcount(aliphatic) == 0

    carbonyl = parse(SMARTS, "[CX3]=[OX1]")
    @test getatom(carbonyl, 1).query == (:and => (
        :and => (:atomsymbol => :C, :isaromatic => false),
        :connectivity => 3
    ))
    @test getatom(carbonyl, 2).query == (:and => (
        :and => (:atomsymbol => :O, :isaromatic => false),
        :connectivity => 1
    ))

    ether = parse(SMARTS, "[#6][OD2][#6]")
    @test getatom(ether, 1).query == (:atomsymbol => :C)
    @test getatom(ether, 2).query == (:and => (
        :and => (:atomsymbol => :O, :isaromatic => false),
        :nodedegree => 2
    ))

    notH = parse(SMARTS, "[!#1]")
    @test getatom(notH, 1).query == (:not => (:atomsymbol => :H))
end

end # smarts.SMARTS
