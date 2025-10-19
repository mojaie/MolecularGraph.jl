#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "isomorphism_clique" begin

@testset "mcs" begin
    # default_logger = global_logger(ConsoleLogger(stdout, Logging.Debug))
    nullg = mcis_constraints(Val(:connection), SimpleGraph())
    @test length(maximum_common_subgraph(nullg, nullg)[1]) == 0
    nullg = mces_constraints(Val(:connection), SimpleGraph())
    @test length(maximum_common_subgraph(nullg, nullg)[1]) == 0

    g = SimpleGraph(Edge.([(1,2), (2,3), (2,4), (2,5), (5,6)]))
    h = SimpleGraph(Edge.([(1,2), (2,3), (3,4), (3,5), (3,6)]))
    gi = mcis_constraints(Val(:connection), g)
    hi = mcis_constraints(Val(:connection), h)
    @test length(maximum_common_subgraph(gi, hi)[1]) == 6
    ge = mces_constraints(Val(:connection), g)
    he = mces_constraints(Val(:connection), h)
    @test length(maximum_common_subgraph(ge, he)[1]) == 5

    g = path_graph(4)
    h = cycle_graph(4)
    gi = mcis_constraints(Val(:connection), g)
    hi = mcis_constraints(Val(:connection), h)
    @test length(maximum_common_subgraph(gi, hi)[1]) == 3
    ge = mces_constraints(Val(:connection), g)
    he = mces_constraints(Val(:connection), h)
    @test length(maximum_common_subgraph(ge, he)[1]) == 3

    g = complete_graph(6)
    h = copy(g)
    rem_edge!(h, Edge(5 => 6))
    gi = mcis_constraints(Val(:connection), g)
    hi = mcis_constraints(Val(:connection), h)
    @test length(maximum_common_subgraph(gi, hi)[1]) == 5
    hij = JSON.parse(JSON.json(hi), ConstraintArrayMCIS{Int,Bool,Int,Int})
    @test length(maximum_common_subgraph(gi, hij)[1]) == 5
    ge = mces_constraints(Val(:connection), g)
    he = mces_constraints(Val(:connection), h)
    @test length(maximum_common_subgraph(ge, he)[1]) == 14
    hej = JSON.parse(JSON.json(he), ConstraintArrayMCES{Int,Bool,Int,Int})
    @test length(maximum_common_subgraph(ge, hej)[1]) == 14

    g = complete_graph(3)
    h = star_graph(4)
    gi = mcis_constraints(Val(:connection), g)
    hi = mcis_constraints(Val(:connection), h)
    @test length(maximum_common_subgraph(gi, hi)[1]) == 2
    ge = mces_constraints(Val(:connection), g)
    he = mces_constraints(Val(:connection), h)
    @test length(maximum_common_subgraph(ge, he)[1]) == 2  # delta-Y transformation
    @test length(maximum_common_subgraph(he, ge)[1]) == 2  # delta-Y transformation
    # global_logger(default_logger)
end

@testset "connected" begin
    # default_logger = global_logger(ConsoleLogger(stdout, Logging.Debug))
    g = disjoint_union(complete_graph(5), complete_graph(4))
    h = disjoint_union(complete_graph(4), complete_graph(3))
    gi = mcis_constraints(Val(:connection), g)
    hi = mcis_constraints(Val(:connection), h)
    @test length(maximum_common_subgraph(gi, hi)[1]) == 7
    @test length(maximum_common_subgraph(gi, hi, connected=true)[1]) == 4
    # global_logger(default_logger)
end

end # isomorphism_clique
