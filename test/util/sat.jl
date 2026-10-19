
@testset "util.sat" begin
    @test is_satisfiable(Vector{Int}[])  # empty formula
    @test !is_satisfiable([Int[]])  # empty clause
    @test is_satisfiable([[1]])
    @test !is_satisfiable([[1], [-1]])
    @test is_satisfiable([[1, 2], [-1, 2], [1, -2]])
    @test !is_satisfiable([[1, 2], [-1, 2], [1, -2], [-1, -2]])
    # unit propagation chain: 1, 1 -> 2, 2 -> 3, !3
    @test !is_satisfiable([[1], [-1, 2], [-2, 3], [-3]])
    # pigeonhole principle: 3 pigeons into 2 holes (unsatisfiable, requires branching)
    p(i, j) = 2 * (i - 1) + j  # pigeon i in hole j
    php = [[p(i, 1), p(i, 2)] for i in 1:3]
    for j in 1:2, i in 1:3, k in (i + 1):3
        push!(php, [-p(i, j), -p(k, j)])
    end
    @test !is_satisfiable(php)
    @test is_satisfiable(php[[1, 2, 4, 5, 6, 7, 8, 9]])  # 2 pigeons
end
