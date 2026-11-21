#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "isomorphism_clique" begin

@testset "constraints" begin
    pairs, arr = connection_constraints(path_graph(3))
    @test pairs == [(1, 2), (1, 3), (2, 3)]
    @test arr == Bool[1, 0, 1]
    pairs, arr = connection_constraints(Int, path_graph(3))
    @test arr == [1, 0, 1]

    pairs, arr = shortest_path_constraints(path_graph(4))
    @test pairs == [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
    @test arr == UInt8[1, 2, 3, 1, 2, 1]
    pairs, arr = shortest_path_constraints(path_graph(4), diameter=2)
    @test pairs == [(1, 2), (1, 3), (2, 3), (2, 4), (3, 4)]
    @test arr == UInt8[1, 2, 1, 2, 1]
    # disconnected pairs are not included
    pairs, arr = shortest_path_constraints(disjoint_union(path_graph(2), path_graph(2)))
    @test pairs == [(1, 2), (3, 4)]

    pairs, arr = any_path_constraints(cycle_graph(4))
    @test pairs == [(1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
    @test arr == UInt32[1, 2, 1, 1, 2, 1]  # bit set of distances (1 if adjacent)
    pairs, arr = any_path_constraints(cycle_graph(6))
    @test arr[findfirst(==((1, 3)), pairs)] == UInt32(0b1010)  # distance 2 and 4
    @test arr[findfirst(==((1, 4)), pairs)] == UInt32(0b100)  # distance 3
    pairs, arr = any_path_constraints(cycle_graph(6), diameter=3)
    @test arr[findfirst(==((1, 3)), pairs)] == UInt32(0b10)
    pairs, arr = any_path_constraints(disjoint_union(path_graph(2), path_graph(2)))
    @test pairs == [(1, 2), (3, 4)]  # disconnected pairs are not included
    @test arr == UInt32[1, 1]

    gi = mcis_constraints(Val(:connection), path_graph(4))
    @test gi isa ConstraintArrayMCIS{Int,Bool,Int,Int}
    @test gi.nv == 4
    @test length(gi.pairs) == length(gi.distances) == length(gi.vattrs1) == length(gi.eattrs) == 6
    gi = mcis_constraints(Val(:shortest), path_graph(4), vmatchvec=v -> 10v, ematchvec=(u, v) -> u + v)
    @test gi isa ConstraintArrayMCIS{Int,UInt8,Int,Int}
    @test gi.vattrs1 == [10 * u for (u, v) in gi.pairs]
    @test gi.vattrs2 == [10 * v for (u, v) in gi.pairs]
    @test gi.eattrs == [u + v for (u, v) in gi.pairs]
    @test mcis_constraints(Val(:any), path_graph(4)) isa ConstraintArrayMCIS{Int,UInt32,Int,Int}

    ge = mces_constraints(Val(:connection), path_graph(4))
    @test ge isa ConstraintArrayMCES{Int,Bool,Int,Int}
    @test ge.nv == 3  # line graph
    @test issetequal(values(ge.revmap), edges(path_graph(4)))
    @test isempty(ge.delta_edges)
    @test isempty(ge.y_edges)
    ge = mces_constraints(Val(:shortest), complete_graph(3))
    @test ge isa ConstraintArrayMCES{Int,UInt8,Int,Int}
    @test length(unique(Set.(ge.delta_edges))) == 1
    ge = mces_constraints(Val(:any), star_graph(4))
    @test ge isa ConstraintArrayMCES{Int,UInt32,Int,Int}
    @test length(ge.y_edges) == 1
end

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

@testset "mcs_constraint_types" begin
    g = SimpleGraph(Edge.([(1,2), (2,3), (2,4), (2,5), (5,6)]))
    h = SimpleGraph(Edge.([(1,2), (2,3), (3,4), (3,5), (3,6)]))
    for c in (:connection, :shortest, :any)
        gi = mcis_constraints(Val(c), g)
        hi = mcis_constraints(Val(c), h)
        mapping, status = maximum_common_subgraph(gi, hi)
        @test status == :done
        @test length(mapping) == 6
        ge = mces_constraints(Val(c), g)
        he = mces_constraints(Val(c), h)
        @test length(maximum_common_subgraph(ge, he)[1]) == 5
    end

    # topological constraint with tolerance
    c6 = mcis_constraints(Val(:shortest), cycle_graph(6))
    p6 = mcis_constraints(Val(:shortest), path_graph(6))
    @test length(maximum_common_subgraph(c6, p6)[1]) == 4
    @test length(maximum_common_subgraph(c6, p6, tolerance=2)[1]) == 5
    c6 = mcis_constraints(Val(:connection), cycle_graph(6))
    p6 = mcis_constraints(Val(:connection), path_graph(6))
    @test length(maximum_common_subgraph(c6, p6)[1]) == 5
    @test_throws ErrorException maximum_common_subgraph(
        mcis_constraints(Val(:shortest), cycle_graph(6)),
        mcis_constraints(Val(:shortest), path_graph(6)), connected=true)
end

@testset "mcs_attributes" begin
    label = [1, 2, 1, 2, 1, 2]
    gi = mcis_constraints(Val(:connection), cycle_graph(6), vmatchvec=v -> label[v])
    hi = mcis_constraints(Val(:connection), path_graph(4), vmatchvec=v -> [1, 2, 1, 2][v])
    mapping, _ = maximum_common_subgraph(gi, hi)
    @test length(mapping) == 4
    @test all(label[k] == [1, 2, 1, 2][v] for (k, v) in mapping)
    hi = mcis_constraints(Val(:connection), path_graph(4), vmatchvec=v -> [2, 2, 1, 1][v])
    mapping, _ = maximum_common_subgraph(gi, hi)
    @test length(mapping) == 2
    @test all(label[k] == [2, 2, 1, 1][v] for (k, v) in mapping)
end

@testset "mcs_options" begin
    k = mcis_constraints(Val(:connection), complete_graph(8))
    mapping, status = maximum_common_subgraph(k, k, targetsize=3)
    @test status == :targetreached
    @test length(mapping) == 3
    k = mces_constraints(Val(:connection), complete_graph(5))
    mapping, status = maximum_common_subgraph(k, k, targetsize=3)
    @test status == :targetreached
    @test length(mapping) <= 3
end

@testset "modular_product_constraints" begin
    gi = mcis_constraints(Val(:connection), path_graph(4))
    @test nv(modular_product(gi, gi)) == 16
    @test ne(modular_product(gi, gi)) == 36
    hi = mcis_constraints(Val(:connection), path_graph(4), vmatchvec=v -> v)
    @test ne(modular_product(hi, hi)) == 6  # only identity mapping
end

@testset "connected" begin
    # default_logger = global_logger(ConsoleLogger(stdout, Logging.Debug))
    g = disjoint_union(complete_graph(5), complete_graph(4))
    h = disjoint_union(complete_graph(4), complete_graph(3))
    gi = mcis_constraints(Val(:connection), g)
    hi = mcis_constraints(Val(:connection), h)
    @test length(maximum_common_subgraph(gi, hi)[1]) == 7
    @test length(maximum_common_subgraph(gi, hi, connected=true)[1]) == 4
    # topological constraints have no pairs between different components
    gi = mcis_constraints(Val(:any), g)
    hi = mcis_constraints(Val(:any), h)
    @test length(maximum_common_subgraph(gi, hi)[1]) == 4

    g = disjoint_union(complete_graph(4), path_graph(4))
    h = disjoint_union(path_graph(3), path_graph(3))
    ge = mces_constraints(Val(:connection), g)
    he = mces_constraints(Val(:connection), h)
    mapping, status = maximum_common_subgraph(ge, he)
    @test status == :done
    @test length(mapping) == 4
    @test issubset(keys(mapping), edges(g))
    @test issubset(values(mapping), edges(h))
    mapping, status = maximum_common_subgraph(ge, he, connected=true)
    @test status == :done
    @test length(mapping) == 2
    # global_logger(default_logger)
end

end # isomorphism_clique
