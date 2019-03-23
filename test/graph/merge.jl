#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "graph.merge" begin

@testset "shallowmerge" begin
    G = pathgraph(5)
    H = pathgraph(5)
    U, nmap, emap = shallowmerge(G, H)
    @test nodecount(U) == 10
    @test edgecount(U) == 8
    @test issetequal(
        [tuple(sort(collect(c))...) for c in connected_components(U)],
        [tuple(1:5...), tuple(6:10...)]
    )
    U, nmap, emap = shallowmerge(G, G)
    @test getnode(U, 4) === getnode(U, 9)
    @test getedge(U, 4) !== getedge(U, 8)
    # TODO: shallowmerge does not copy Node and Edge objects
end

end # graph.merge
