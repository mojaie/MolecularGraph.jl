#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "isomorphism_clique" begin

@testset "mcs" begin
    # default_logger = global_logger(ConsoleLogger(stdout, Logging.Debug))
    nullg = SimpleGraph()
    @test size(maximum_common_subgraph(nullg, nullg)) == 0
    @test size(maximum_common_edge_subgraph(nullg, nullg)) == 0

    g = SimpleGraph(Edge.([(1,2), (2,3), (2,4), (2,5), (5,6)]))
    h = SimpleGraph(Edge.([(1,2), (2,3), (3,4), (3,5), (3,6)]))
    @test size(maximum_common_subgraph(g, h)) == 6
    @test size(maximum_common_edge_subgraph(g, h)) == 5

    g = path_graph(4)
    h = cycle_graph(4)
    @test size(maximum_common_subgraph(g, h)) == 3
    @test size(maximum_common_edge_subgraph(g, h)) == 3

    g = complete_graph(6)
    h = copy(g)
    rem_edge!(h, Edge(5 => 6))
    @test size(maximum_common_subgraph(g, h)) == 5
    @test size(maximum_common_edge_subgraph(g, h)) == 14

    g = complete_graph(3)
    h = star_graph(4)
    @test size(maximum_common_subgraph(g, h)) == 2
    @test size(maximum_common_edge_subgraph(g, h)) == 2  # delta-Y transformation
    @test size(maximum_common_edge_subgraph(h, g)) == 2  # delta-Y transformation
    # global_logger(default_logger)
end

@testset "connected" begin
    # default_logger = global_logger(ConsoleLogger(stdout, Logging.Debug))
    g, vmap = disjoint_union(complete_graph(5), complete_graph(4))
    h, vmap = disjoint_union(complete_graph(4), complete_graph(3))
    @test size(maximum_common_subgraph(g, h)) == 7
    @test size(maximum_common_subgraph(g, h, connected=true)) == 4
    # global_logger(default_logger)
end

end # isomorphism_clique
