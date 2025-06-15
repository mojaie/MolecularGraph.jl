#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "isomorphism_clique" begin

@testset "mcs" begin
    # default_logger = global_logger(ConsoleLogger(stdout, Logging.Debug))
    nullg = SimpleGraph()
    @test length(maximum_common_subgraph(nullg, nullg)[1]) == 0
    @test length(maximum_common_edge_subgraph(nullg, nullg)[1]) == 0

    g = SimpleGraph(Edge.([(1,2), (2,3), (2,4), (2,5), (5,6)]))
    h = SimpleGraph(Edge.([(1,2), (2,3), (3,4), (3,5), (3,6)]))
    @test length(maximum_common_subgraph(g, h)[1]) == 6
    @test length(maximum_common_edge_subgraph(g, h)[1]) == 5

    g = path_graph(4)
    h = cycle_graph(4)
    @test length(maximum_common_subgraph(g, h)[1]) == 3
    @test length(maximum_common_edge_subgraph(g, h)[1]) == 3

    g = complete_graph(6)
    h = copy(g)
    rem_edge!(h, Edge(5 => 6))
    @test length(maximum_common_subgraph(g, h)[1]) == 5
    @test length(maximum_common_edge_subgraph(g, h)[1]) == 14

    g = complete_graph(3)
    h = star_graph(4)
    @test length(maximum_common_subgraph(g, h)[1]) == 2
    @test length(maximum_common_edge_subgraph(g, h)[1]) == 2  # delta-Y transformation
    @test length(maximum_common_edge_subgraph(h, g)[1]) == 2  # delta-Y transformation
    # global_logger(default_logger)
end

@testset "connected" begin
    # default_logger = global_logger(ConsoleLogger(stdout, Logging.Debug))
    g = disjoint_union(complete_graph(5), complete_graph(4))
    h = disjoint_union(complete_graph(4), complete_graph(3))
    @test length(maximum_common_subgraph(g, h)[1]) == 7
    @test length(maximum_common_subgraph(g, h, connected=true)[1]) == 4
    # global_logger(default_logger)
end

end # isomorphism_clique
