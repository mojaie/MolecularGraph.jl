#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "clique" begin

@testset "max_clique" begin
    nullg = SimpleGraph()
    @test maximum_clique(nullg)[1] == []

    noedges = SimpleGraph(5)
    @test length(all_maximal_cliques(noedges)[1]) == 5
    @test length(maximum_clique(noedges)[1]) == 1

    p7 = path_graph(7)
    @test length(all_maximal_cliques(p7)[1]) == 6
    @test length(maximum_clique(p7)[1]) == 2

    k5 = complete_graph(5)
    @test length(all_maximal_cliques(k5)[1]) == 1
    @test length(maximum_clique(k5)[1]) == 5

    w8 = wheel_graph(8)
    @test length(all_maximal_cliques(w8)[1]) == 7
    @test length(maximum_clique(w8)[1]) == 3

    petersen = smallgraph(:petersen)
    @test length(all_maximal_cliques(petersen)[1]) == 15
    @test length(maximum_clique(petersen)[1]) == 2

    karate = smallgraph(:karate)
    @test length(all_maximal_cliques(karate)[1]) == 36
    @test length(maximum_clique(karate)[1]) == 5

end

@testset "max_conn_clique" begin
    nullg = SimpleGraph()
    @test maximum_conn_clique(nullg)[1] == []

    noedges = SimpleGraph(5)
    @test length(all_maximal_conn_cliques(noedges)[1]) == 5
    @test length(maximum_conn_clique(noedges)[1]) == 1

    k5 = complete_graph(5)
    @test length(all_maximal_conn_cliques(k5)[1]) == 5
    @test length(maximum_conn_clique(k5)[1]) == 1
    connfunc = n -> Set.(Dict(
        1 => ([2], [3, 4, 5]), 2 => ([1, 3], [4, 5]), 3 => ([2], [1, 4, 5]),
        4 => ([5], [1, 2, 3]), 5 => ([4], [1, 2, 3]))[n])
    @test length(all_maximal_conn_cliques(k5, connfunc=connfunc)[1]) == 2
    @test length(maximum_conn_clique(k5, connfunc=connfunc)[1]) == 3
    connfunc = n -> (Set(neighbors(k5, n)), Set{Int}())
    @test length(all_maximal_conn_cliques(k5, connfunc=connfunc)[1]) == 1
    @test length(maximum_conn_clique(k5, connfunc=connfunc)[1]) == 5
end

end # clique