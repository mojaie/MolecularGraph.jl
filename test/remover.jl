
@testset "remover" begin

@testset "hydrogen" begin
    ethanol = parsesmiles("[H]C([H])([H])C([H])([H])O", mutable=true)
    @test atomcount(ethanol) == 8
    remove_H!(ethanol)
    @test atomcount(ethanol) == 3
end

@testset "water" begin
    CuSO4_5H2O = parsesmiles("[Cu2+].[O-]S(=O)(=O)[O-].O.O.O.O.O", mutable=true)
    remove_H!(CuSO4_5H2O)
    @test atomcount(CuSO4_5H2O) == 11
    CuSO4 = removewater(CuSO4_5H2O)
    @test atomcount(CuSO4) == 6
end

end # substructure
