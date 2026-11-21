#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "bipartite" begin

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
    @test maxcard(Set(1:3), Set(1:3), Dict(i => Set{Int}() for i in 1:3)) == 0
end

@testset "maxcardmap" begin
    relation = Dict(1 => Set([1]), 2 => Set([1, 2]), 3 => Set([2, 3]))
    @test maxcardmap(Set(1:3), Set(1:3), relation) == Dict(1 => 1, 2 => 2, 3 => 3)
    relation = Dict(
        1 => Set([1, 2]), 2 => Set([3]),
        3 => Set([4, 5]), 4 => Set([3, 4]), 5 => Set([2, 5])
    )
    mapping = maxcardmap(Set(1:5), Set(1:5), relation)
    @test length(mapping) == 5
    @test all(j in relation[i] for (i, j) in mapping)  # valid adjacency
    @test allunique(values(mapping))
    relation = Dict(1 => Set(1:6), 2 => Set([4]), 3 => Set([4]), 4 => Set([4]))
    mapping = maxcardmap(Set(1:4), Set(1:6), relation)
    @test length(mapping) == 2
    @test all(j in relation[i] for (i, j) in mapping)
    @test allunique(values(mapping))
    # adjacency given by a function
    @test maxcardmap(1:3, 1:3, (i, j) -> i == j) == Dict(1 => 1, 2 => 2, 3 => 3)
end

end # bipartite
