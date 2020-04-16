
@testset "util.math" begin

@testset "factorial" begin
    @test_throws DomainError logfactorial(-1)
    @test isapprox(logfactorial(0), 0.0, atol=1e-3)
    @test isapprox(logfactorial(1), 0.0, atol=1e-3)
    @test isapprox(logfactorial(100), 363.739, atol=1e-3)
end

end # util.math
