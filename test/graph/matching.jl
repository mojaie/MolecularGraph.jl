#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "matching" begin

@testset "max_matching" begin
    # default_logger = global_logger(ConsoleLogger(stdout, Logging.Debug))
    p7 = path_graph(7)
    @test !is_perfect_matching(p7)
    p8 = path_graph(8)
    @test is_perfect_matching(p8)
    c5 = cycle_graph(5)
    @test !is_perfect_matching(c5)
    c6 = cycle_graph(6)
    @test is_perfect_matching(c6)
    s10 = star_graph(10)
    @test length(max_matching(s10)) == 1
    k7 = complete_graph(7)
    @test length(max_matching(k7)) == 3
    k47 = complete_bipartite_graph(4, 7)
    @test length(max_matching(k47)) == 4
    l8 = ladder_graph(8)
    @test is_perfect_matching(l8)
    dodeca = smallgraph(:dodecahedral)
    @test length(max_matching(dodeca)) == 10
    k20 = complete_graph(20)
    @test length(max_matching(k20)) == 10
    # global_logger(default_logger)
end

end # matching
