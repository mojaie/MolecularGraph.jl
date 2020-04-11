
@testset "mass" begin

@testset "mass" begin
    iron = SDFileAtom(:Fe)
    @test atomnumber(iron) == 26
    @test standardweight(iron) == (55.845, 0.002)
    @test standardweight(Float64, iron) == 55.84
    @test monoisotopicmass(iron) == (55.93493633, 4.9e-7)
    @test monoisotopicmass(Float64, iron) == 55.934936
    @test exactmass(iron) == (55.93493633, 4.9e-7)
    @test exactmass(Float64, iron) == 55.934936
    
    og = SDFileAtom(:Og)
    @test atomnumber(og) == 118
    @test standardweight(og) === (NaN, NaN)
    @test exactmass(og) === (NaN, NaN)
    og294 = SDFileAtom(:Og, 0, 1, 294, [0.0, 0.0, 0.0])
    @test monoisotopicmass(og294) === (NaN, NaN)
    @test exactmass(og294) === (294.21392, 0.00071)

    etoh = smilestomol("CCO")
    wt = standardweight(etoh)
    @test isapprox(wt[1], 46.069, atol=1e-3)
    @test isapprox(wt[2], 0.004, atol=1e-3)
    @test standardweight(Float64, etoh) == 46.07
    iso = monoisotopicmass(etoh)
    @test isapprox(iso[1], 46.0418648130, atol=1e-10)
    @test isapprox(iso[2], 7.1e-10, atol=1e-10)
    ex = exactmass(etoh)
    @test isapprox(iso[1], 46.0418648130, atol=1e-10)
    @test isapprox(iso[2], 7.1e-10, atol=1e-10)
end

end  # mass
