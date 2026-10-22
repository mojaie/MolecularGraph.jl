
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

@testset "deuterium" begin
    @test nominal_mass(SMILESAtom(:D)) == 2
    @test nominal_mass(SMILESAtom(:T)) == 3
    @test isapprox(monoiso_mass(SMILESAtom(:D)), 2.01410177812, atol=1e-10)
    @test isapprox(standard_weight(SMILESAtom(:D)), 2.01410177812, atol=1e-10)
    @test isapprox(standard_weight(SMILESAtom(:T)), 3.0160492779, atol=1e-10)
    @test isapprox(standard_weight(SMILESAtom(:D, 0, 0, 2)), 2.01410177812, atol=1e-10)
    @test isapprox(standard_weight(smilestomol("[2H]O[2H]")), standard_weight(smilestomol("O")) + 2 * (2.01410177812 - 1.008), atol=1e-10)
    @test MolecularGraph.atom_markup(:D, 2, 0) == [[(:default, "D")]]
    @test MolecularGraph.atom_markup(:D, 0, 0) == [[(:default, "D")]]
    @test MolecularGraph.atom_markup(:T, 3, 0) == [[(:default, "T")]]
    @test MolecularGraph.atom_markup(:T, 0, 0) == [[(:default, "T")]]
    @test MolecularGraph.atom_markup(:H, 2, 0) == [[(:sup, "2"), (:default, "H")]]
end

@testset "isotopic_composition" begin
    @test length(isotopic_composition(:C, 100; threshold=0.01)) == 5
    @test length(isotopic_composition(:H, 1000; threshold=0.01)) == 2
    @test length(isotopic_composition(smilestomol("CCl"))) == 4
    # sorted by mass in ascending order
    @test issorted(isotopic_composition(:C, 100; threshold=0.01), by=x->x[1])
    @test issorted(isotopic_composition(smilestomol("CCl")), by=x->x[1])
    # display(simulatemassspec(smilestomol("CCl")))
    # keyword arguments are passed to the corresponding functions
    spec = simulate_massspec(smilestomol("CCl"))
    @test size(spec, 2) == 2
    @test size(simulate_massspec(smilestomol("CCl"); threshold=0.01, rate=0.01), 1) < size(spec, 1)
end

end  # mass
