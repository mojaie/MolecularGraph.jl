
using MolecularGraph: lookahead, forward!, backtrack!


@testset "smarts.base" begin

@testset "base" begin
    state1 = SmilesParser("C1SC(C=O)CCC1", true)
    state2 = SmartsParser(raw"*OC$([Cv4H2+0])", false)
    @test read(state1) == 'C'
    @test lookahead(state1, 2) == 'S'
    @test read(state2) == '*'
    forward!(state2)
    @test read(state2) == 'O'
    forward!(state2, 2)
    @test read(state2) == '$'
    backtrack!(state2)
    @test lookahead(state2, 1) == '$'
    @test !state2.done
    forward!(state2, 13)
    @test state2.done
end

@testset "equivalence" begin
    q = (:any => true)
    @test isequivalent(q, q)
    @test query_contains(q, q)

    q1 = (:and => (:atomsymbol => :C, :isaromatic => false, :mass => 14))
    q2 = (:and => (:mass => 14, :atomsymbol => :C, :isaromatic => false))
    @test isequivalent(q1, q2)
    @test query_contains(q1, q2)
    q3 = (:and => (:atomsymbol => :C, :isaromatic => false, :mass => 13))
    @test !isequivalent(q1, q3)
    @test !query_contains(q1, q3)
    q4 = (:and => (:atomsymbol => :C, :isaromatic => false))
    @test !isequivalent(q1, q4)
    @test !query_contains(q1, q4)
    @test query_contains(q4, q1)
    q5 = (:and => (:atomsymbol => :C, :hydrogenconnected => 2, :isaromatic => false, :mass => 14))
    @test !isequivalent(q1, q5)
    @test query_contains(q1, q5)
    @test !query_contains(q5, q1)

    q6 = (:or => (:atomsymbol => :C, :isaromatic => false, :mass => 14))
    @test !isequivalent(q1, q6)
    @test !query_contains(q1, q6)
    q7 = (:or => (:atomsymbol => :C, :isaromatic => false))
    @test !isequivalent(q6, q7)
    @test query_contains(q6, q7)
    @test !query_contains(q7, q6)
    q8 = (:or => (:atomsymbol => :C, :hydrogenconnected => 2, :isaromatic => false, :mass => 14))
    @test !isequivalent(q7, q8)
    @test !query_contains(q7, q8)
    @test query_contains(q8, q7)
    q9 = :atomsymbol => :C
    @test !isequivalent(q8, q9)
    @test query_contains(q8, q9)
    @test !query_contains(q9, q8)

    nested1 = (:and => (
        :or => (:atomsymbol => :O, :atomsymbol => :S),
        :or => (:hydrogenconnected => 1, :hydrogenconnected => 2),
        :isaromatic => false
    ))
    nested2 = (:and => (
        :or => (:atomsymbol => :O, :atomsymbol => :S),
        :hydrogenconnected => 2,
        :isaromatic => false
    ))
    nested3 = (:and => (
        :or => (:atomsymbol => :O, :atomsymbol => :S),
        :or => (:hydrogenconnected => 1, :not => (:hydrogenconnected => 2)),
        :isaromatic => false
    ))
    @test !isequivalent(nested1, nested2)
    @test query_contains(nested1, nested2)
    @test !query_contains(nested2, nested1)
    @test !isequivalent(nested1, nested3)
    @test !query_contains(nested1, nested3)

    nested4 = (:or => (
        :and => (:atomsymbol => :S, :connectivity => 3, :charge => 0),
        :and => (:atomsymbol => :S, :connectivity => 3, :charge => 1),
        :and => (:atomsymbol => :O, :isaromatic => false)
    ))
    nested5 = (:or => (
        :and => (:atomsymbol => :S, :charge => 0),
        :and => (:atomsymbol => :S, :charge => 1),
        :and => (:atomsymbol => :O, :isaromatic => false)
    ))
    nested6 = (:or => (
        :and => (:atomsymbol => :S, :charge => 0),
        :and => (:atomsymbol => :S, :charge => 1),
        :and => (:atomsymbol => :O, :isaromatic => true)
    ))
    nested7 = (:or => (
        :and => (:atomsymbol => :S, :connectivity => 3, :charge => 0),
        :and => (:atomsymbol => :S, :charge => 1, :connectivity => 3),
        :and => (:atomsymbol => :O, :isaromatic => false)
    ))
    @test !isequivalent(nested4, nested5)
    @test !query_contains(nested4, nested5)
    @test query_contains(nested5, nested4)
    @test !isequivalent(nested4, nested6)
    @test !query_contains(nested6, nested4)
    @test isequivalent(nested4, nested7)
end

end # smarts.base
