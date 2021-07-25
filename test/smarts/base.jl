
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

end # smarts.base
