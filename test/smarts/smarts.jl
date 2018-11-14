
@testset "smarts.smarts" begin

@testset "chain" begin
    nullmol = parse(SMARTS, "")
    @test atomcount(nullmol) == 0
    @test bondcount(nullmol) == 0

    aliphatic = parse(SMARTS, "C")
    @test getatom(aliphatic, 1).query == (
        :and => (:Symbol => :C, :Aromatic => false)
    )
    @test bondcount(aliphatic) == 0

    carbonyl = parse(SMARTS, "[CX3]=[OX1]")
    @test getatom(carbonyl, 1).query == (:and => (
        :and => (:Symbol => :C, :Aromatic => false),
        :Connectivity => 3
    ))
    @test getatom(carbonyl, 2).query == (:and => (
        :and => (:Symbol => :O, :Aromatic => false),
        :Connectivity => 1
    ))

    ether = parse(SMARTS, "[#6][OD2][#6]")
    @test getatom(ether, 1).query == (:Symbol => :C)
    @test getatom(ether, 2).query == (:and => (
        :and => (:Symbol => :O, :Aromatic => false),
        :Degree => 2
    ))

    notH = parse(SMARTS, "[!#1]")
    @test getatom(notH, 1).query == (:not => (:Symbol => :H))
end

end # smarts.SMARTS
