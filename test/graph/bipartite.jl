#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.bipartite" begin

@testset "maxcardmatch" begin
    relation1 = Dict(
        1 => Set([1]), 2 => Set([1, 2, 3]),
        3 => Set([1, 2]), 4 => Set([2, 3, 4])
    )
    @test maxcard(Set(1:4), Set(1:4), relation1) == 4
    relation2 = Dict(
        1 => Set([2]), 2 => Set([2, 3]),
        3 => Set(1:5), 4 => Set([3]), 5 => Set([4])
    )
    @test maxcard(Set(1:5), Set(1:5), relation2) == 4
    relation3 = Dict(
        1 => Set([1, 2]), 2 => Set([3]),
        3 => Set([4, 5]), 4 => Set([3, 4]), 5 => Set([2, 5])
    )
    @test maxcard(Set(1:5), Set(1:5), relation3) == 5
    relation4 = Dict(
        1 => Set(1:6), 2 => Set([4]), 3 => Set([4]), 4 => Set([4])
    )
    @test maxcard(Set(1:4), Set(1:6), relation4) == 2
end

@testset "twocoloring" begin
    K5_5 = completebipartite(5, 5)
    a, b = twocoloring(K5_5)
    @test length(a) == 5
    @test length(b) == 5

    CL5 = circularladder(5)
    n = twocoloring(CL5)
    @test n === nothing

    CL6 = circularladder(6)
    a, b = twocoloring(CL6)
    @test length(a) == 6
    @test length(b) == 6
end

end # graph.bipartite
