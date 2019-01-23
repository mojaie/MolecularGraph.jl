
@testset "util.iterator" begin

@testset "combinations" begin
    @test length(collect(combinations(1:10, 1))) == 10
    @test length(collect(combinations(1:10, 2))) == 45
    @test length(collect(combinations(1:10, 3))) == 120
    @test length(collect(combinations(1:10, 8))) == 45
    @test length(collect(combinations(1:10, 9))) == 10
    @test length(collect(combinations(1:10, 10))) == 1
    @test length(collect(combinations(1:10, 0))) == 1
    @test_throws DomainError combinations(1:10, -1)
    @test_throws ErrorException combinations(1:10, 11)
    @test_throws ErrorException combinations([], 1)
    @test_throws MethodError combinations(1:10, 1.5)
end

end # util.iterator
