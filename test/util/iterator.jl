
@testset "util.iterator" begin

@testset "combinations" begin
    @test length(collect(combinations(10, 1))) == 10
    @test length(collect(combinations(10, 2))) == 45
    @test length(collect(combinations(10, 3))) == 120
    @test length(collect(combinations(10, 8))) == 45
    @test length(collect(combinations(10, 9))) == 10
    @test length(collect(combinations(10, 10))) == 1
    @test length(collect(combinations(10, 0))) == 1
    @test_throws DomainError combinations(10, -1)
    @test_throws ErrorException combinations(10, 11)
    @test_throws ErrorException combinations(0, 1)
    @test_throws MethodError combinations(10, 1.5)
end

@testset "sortstablemax" begin
    @test sortstablemax([1, 4, 2, 7, 5]) == 7
    @test sortstablemax([8:10, 1:3, [], 4:23, 1:9], by=length) == 4:23
    @test sortstablemax([8:10, 1:3, 3:22, 4:23, 1:9], by=length) == 3:22
    @test_throws ErrorException sortstablemax([])
    @test sortstablemax([], init=0) == 0
end

@testset "sortstablemin" begin
    @test sortstablemin([1, 4, 2, 7, 5]) == 1
    @test sortstablemin([8:10, 1:3, [], 4:23, 1:9], by=length) == []
    @test sortstablemin([8:10, 1:3, 3:22, 4:23, 1:9], by=length) == 8:10
    @test_throws ErrorException sortstablemin([])
    @test sortstablemin([], init=0) == 0
end
end # util.iterator
