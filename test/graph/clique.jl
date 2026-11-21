#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "clique" begin

@testset "max_clique" begin
    nullg = SimpleGraph()
    @test maximum_clique_mg(nullg)[1] == []

    noedges = SimpleGraph(5)
    @test length(all_maximal_cliques(noedges)[1]) == 5
    @test length(maximum_clique_mg(noedges)[1]) == 1

    p7 = path_graph(7)
    @test length(all_maximal_cliques(p7)[1]) == 6
    @test length(maximum_clique_mg(p7)[1]) == 2

    k5 = complete_graph(5)
    @test length(all_maximal_cliques(k5)[1]) == 1
    @test length(maximum_clique_mg(k5)[1]) == 5

    w8 = wheel_graph(8)
    @test length(all_maximal_cliques(w8)[1]) == 7
    @test length(maximum_clique_mg(w8)[1]) == 3

    petersen = smallgraph(:petersen)
    @test length(all_maximal_cliques(petersen)[1]) == 15
    @test length(maximum_clique_mg(petersen)[1]) == 2

    karate = smallgraph(:karate)
    @test length(all_maximal_cliques(karate)[1]) == 36
    @test length(maximum_clique_mg(karate)[1]) == 5

    w5 = wheel_graph(5)
    cliques, status = all_maximal_cliques(w5)
    @test status == :done
    @test issetequal(Set.(cliques), Set.([[1, 2, 3], [1, 3, 4], [1, 4, 5], [1, 2, 5]]))
    maxclq, status = maximum_clique_mg(karate)
    @test status == :done
    @test all(has_edge(karate, u, v) for u in maxclq for v in maxclq if u < v)
end

@testset "clique_options" begin
    k6 = complete_graph(6)
    cliques, status = all_maximal_cliques(k6, targetsize=3)
    @test status == :targetreached
    @test length(only(cliques)) == 3
    maxclq, status = maximum_clique_mg(k6, targetsize=3)
    @test status == :targetreached
    @test length(maxclq) == 3
    karate = smallgraph(:karate)
    @test all_maximal_cliques(karate, timeout=0)[2] == :timedout
    @test maximum_clique_mg(karate, timeout=0)[2] == :timedout
    connfunc = n -> (Set(neighbors(karate, n)), Set{Int}())
    @test maximum_conn_clique(karate, connfunc=connfunc, timeout=0)[2] == :timedout
    # custom postprocess
    sizes = Int[]
    all_maximal_cliques(path_graph(4), postprocess=state -> push!(sizes, length(state.Q)))
    @test sizes == [2, 2, 2]
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
    @test length(maximum_conn_clique(k5, connfunc=connfunc, targetsize=3)[1]) == 3
    @test maximum_conn_clique(k5, connfunc=connfunc, targetsize=3)[2] == :targetreached
end

@testset "connfuncgen" begin
    # modular product of P3 and P3: the clique of connected node pairs is the identity mapping
    prod, isconn = modular_product(path_graph(3), path_graph(3))
    connfunc = connfuncgen(prod, isconn)
    @test connfunc(1) == (Set([5]), Set([9]))
    for n in vertices(prod)
        conn, disconn = connfunc(n)
        @test isempty(intersect(conn, disconn))
        @test union(conn, disconn) == Set(neighbors(prod, n))
    end
    maxclq, status = maximum_conn_clique(prod, connfunc=connfunc)
    @test status == :done
    @test length(maxclq) == 3
    # disconnected clique (mapping 1 => 3, 2 => 2, 3 => 1) can also be found
    @test length(maximum_clique_mg(prod)[1]) == 3
    prod, isconn = modular_product(path_graph(3), path_graph(2))
    @test length(maximum_conn_clique(prod, connfunc=connfuncgen(prod, isconn))[1]) == 2
end

@testset "approx_maximum_clique" begin
    @test isempty(approx_maximum_clique(SimpleGraph()))
    @test issetequal(approx_maximum_clique(complete_graph(5)), 1:5)
    @test length(approx_maximum_clique(path_graph(4))) == 2
    @test length(approx_maximum_clique(SimpleGraph(3))) == 1
    for g in (smallgraph(:karate), smallgraph(:petersen), wheel_graph(8))
        clq = approx_maximum_clique(g)
        @test !isempty(clq)
        @test allunique(clq)
        @test all(has_edge(g, u, v) for u in clq for v in clq if u < v)
        @test length(clq) <= length(maximum_clique_mg(g)[1])
    end
end

end # clique