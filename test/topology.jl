
@testset "topology" begin

@testset "resolve_inclusion" begin
    R = resolve_inclusion
    C = canonicalize_cycle
    a = [1, 2, 3, 4, 5]
    b = [1, 2, 3, 4, 5, 6, 7, 8, 9]
    @test C(R(a, b)[2]) == [1, 5, 6, 7, 8, 9]
    b = [1, 6, 7, 8, 9]
    @test R(a, b) == nothing
    b = [3, 4, 5, 6, 7, 8, 9]
    @test R(a, b) == nothing
    b = [2, 3, 4, 5, 6, 7, 8]
    @test C(R(a, b)[2]) == [1, 2, 8, 7, 6, 5]
    a = [3, 4, 5, 6, 7, 8, 9, 10]
    b = [1, 2, 8, 7, 6, 5, 4]
    @test C(R(a, b)[1]) == [1, 2, 8, 9, 10, 3, 4]
    a = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11]
    b = [8, 9, 10, 11, 12]
    @test C(R(a, b)[1]) == [1, 2, 3, 4, 5, 6, 7, 8, 12, 11]
end

"""
@testset "topology" begin
    RESOURCE_DIR = joinpath(dirname(@__FILE__), "..", "_resources", "DrugBank")
    phe = loadsdfmol(open(joinpath(RESOURCE_DIR, "Phe.mol")))
    @test phe.rings[1].arr == UInt16[6, 7, 10, 12, 11, 8]
    @test phe.scaffolds == [[1]]
    @test length(phe.isolated) == 0
    prem = loadsdfmol(open(joinpath(RESOURCE_DIR, "Premarin.mol")))
    @test prem.scaffolds == [[1, 2, 3, 4]]
    @test prem.isolated == [[28]]
end
"""

end # topology
