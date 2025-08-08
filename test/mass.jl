
@testset "mass" begin

@testset "mass" begin
    iron = SDFAtom(;symbol=:Fe)
    @test standard_weight_unc(iron) == (55.845, 0.002)
    @test standard_weight(iron, 2) == 55.84
    @test monoiso_mass_unc(iron) == (55.93493633, 4.9e-7)
    @test monoiso_mass(iron, 6) == 55.934936
    @test exact_mass_unc(iron) == (55.93493633, 4.9e-7)
    @test exact_mass(iron, 6) == 55.934936

    og = SDFAtom(;symbol=:Og)
    @test standard_weight_unc(og) === (NaN, NaN)
    @test exact_mass_unc(og) === (NaN, NaN)
    og294 = SDFAtom(;symbol=:Og, isotope=294)
    @test monoiso_mass_unc(og294) === (NaN, NaN)
    @test exact_mass_unc(og294) === (294.21392, 0.00071)

    etoh = smilestomol("CCO")
    wt, wtunc = standard_weight_unc(etoh)
    @test isapprox(wt, 46.069, atol=1e-3)
    @test isapprox(wtunc, 0.004, atol=1e-3)
    @test standard_weight(etoh, 2) == 46.07
    mo, mounc = monoiso_mass_unc(etoh)
    @test isapprox(mo, 46.0418648130, atol=1e-10)
    @test isapprox(mounc, 7.1e-10, atol=1e-10)
    ex, exunc = exact_mass_unc(etoh)
    @test isapprox(ex, 46.0418648130, atol=1e-10)
    @test isapprox(exunc, 7.1e-10, atol=1e-10)

    etohd6 = smilestomol("[2H]C([2H])([2H])C([2H])([2H])O[2H]")
    @test standard_weight(etohd6, 2) == 52.11
    @test monoiso_mass(etohd6, 6) == 46.041865
    @test exact_mass(etohd6, 6) == 52.079525
end

@testset "isotopic_composition" begin
    @test length(isotopic_composition(:C, 100; threshold=0.01)) == 5
    @test length(isotopic_composition(:H, 1000; threshold=0.01)) == 2
    @test length(isotopic_composition(smilestomol("CCl"))) == 4
    # display(simulatemassspec(smilestomol("CCl")))
end

end  # mass
