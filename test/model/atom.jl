
@testset "model.atom" begin

@testset "atom" begin
    @test atomnumber(:Fe) == 26

    sdfa = SDFileAtom(:Fe, 2, 1, nothing, [1.0, 2.0, 0.0], :anticlockwise)
    @test sdfa.symbol === :Fe
    @test atomnumber(sdfa) == 26
    @test sdfa.charge == 2
    @test sdfa.multiplicity == 1
    @test sdfa.coords == [1.0, 2.0, 0.0]
    @test sdfa.stereo === :anticlockwise

    smia = SmilesAtom(:C, 0, 1, 13.0, false, :clockwise)
    @test smia.symbol === :C
    @test atomnumber(smia) == 6
    @test smia.charge === 0
    @test smia.multiplicity == 1
    @test !smia.isaromatic
    @test smia.stereo === :clockwise
end

@testset "associate_operations" begin
    and1 = (:and => (
        :v => :a,
        :and => (:v => :b, :v => :c)
    ))
    @test Set(associate_operations(and1).second) == Set((:v => :a, :v => :b, :v => :c))

    and2 = (:and => (
        :v => :a,
        :v => :b,
        :and => (:v => :c, :v => :d, :v => :e),
        :or => (:v => :f, :v => :g),
        :not => :v => :h
    ))
    @test Set(associate_operations(and2).second) == Set((
        :v => :a, :v => :b, :v => :c, :v => :d, :v => :e,
        :or => (:v => :f, :v => :g),
        :not => :v => :h
    ))

    or1 = (:or => (
        :or => (:v => :a, :v => :b),
        :v => :c
    ))
    @test Set(associate_operations(or1).second) == Set((:v => :a, :v => :b, :v => :c))

    or2 = (:or => (
        :and => (:v => :a, :v => :b),
        :and => (:v => :c, :v => :d),
        :or => (:v => :e, :v => :f)
    ))
    @test Set(associate_operations(or2).second) == Set((
        :and => (:v => :a, :v => :b),
        :and => (:v => :c, :v => :d),
        :v => :e, :v => :f
    ))

    not1 = (:not => :v => :a)
    @test associate_operations(not1) == (:not => :v => :a)

    nested1 = (:and => (
        :and => (
            :and => (:v => :a, :and => (:v => :b, :v => :c)),
            :v => :d
        ),
        :v => :e
    ))
    @test Set(associate_operations(nested1).second) == Set((:v => :a, :v => :b, :v => :c, :v => :d, :v => :e))
end


@testset "isequivalent" begin
    fml1 = :v => :a
    fml2 = :v => :a
    @test isequivalent(fml1, fml2)
    fml1 = :v => :b
    fml2 = :v => :c
    @test !isequivalent(fml1, fml2)
    fml1 = :x => :b
    fml2 = :v => :b
    @test !isequivalent(fml1, fml2)

    fml1 = :not => :v => :a
    fml2 = :not => :v => :a
    @test isequivalent(fml1, fml2)
    fml1 = :not => :v => :a
    fml2 = :not => :v => :b
    @test !isequivalent(fml1, fml2)

    fml1 = :and => (:v => :a, :v => :b, :v => :c)
    fml2 = :and => (:v => :c, :v => :b, :v => :a)
    @test isequivalent(fml1, fml2)
    fml1 = :or => (:v => :a, :v => :b, :v => :c)
    fml2 = :or => (:v => :c, :v => :b, :v => :a)
    @test isequivalent(fml1, fml2)
    fml1 = :and => (:v => :a, :v => :b, :v => :c)
    fml2 = :and => (:v => :a, :v => :b, :v => :c, :v => :d)
    @test !isequivalent(fml1, fml2)
    fml1 = :or => (:v => :a, :v => :b, :v => :c)
    fml2 = :or => (:x => :a, :y => :b, :z => :c)
    @test !isequivalent(fml1, fml2)

    fml1 = :and => (:v => :a, :or => (:v => :b, :v => :c))
    fml2 = :and => (:v => :a, :or => (:v => :c, :v => :b))
    @test isequivalent(fml1, fml2)
    fml1 = :and => (:v => :a, :and => (:v => :b, :v => :c))
    fml2 = :and => (:v => :a, :or => (:v => :b, :v => :c))
    @test !isequivalent(fml1, fml2)
end

end # atom