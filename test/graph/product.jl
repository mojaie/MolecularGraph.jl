#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

@testset "product" begin

@testset "modularproduct" begin
    g = path_graph(3)
    prod, isconn = modular_product(g, g)
    @test nv(prod) == 9
    @test ne(prod) == 10
end


end # graph.product
