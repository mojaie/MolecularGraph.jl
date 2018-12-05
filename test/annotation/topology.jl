
@testset "annotation.topology" begin

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
