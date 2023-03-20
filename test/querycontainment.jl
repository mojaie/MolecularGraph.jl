
using MolecularGraph:
    resolve_disjoint_not, resolve_recursive, generate_truthtable, querymatch, querypropmap

@testset "querycontainment" begin

@testset "resolve_disjoint_not" begin
    oors = QueryOperator(:or, [
        QueryLiteral(:symbol, :O), QueryLiteral(:symbol, :S)
    ])  # [O,S]
    notn = QueryOperator(:not, [QueryLiteral(:symbol, :N)])  # [!#7]
    oors_ = resolve_disjoint_not(oors, querypropmap(notn))  # [O,S]
    @test oors_.key === :or
    @test oors_.value[1] == QueryLiteral(:symbol, :O)
    @test oors_.value[2] == QueryLiteral(:symbol, :S)
    notn_ = resolve_disjoint_not(notn, querypropmap(oors))  # [!#7,O,S]
    @test notn_.key === :or
    @test notn_.value[1].key == :not
    @test notn_.value[1].value[1] == QueryLiteral(:symbol, :N)
    @test notn_.value[2] == QueryLiteral(:symbol, :O)
    @test notn_.value[3] == QueryLiteral(:symbol, :S)
    noto = QueryOperator(:not, [QueryLiteral(:symbol, :O)])  # [!#6]
    noto_ = resolve_disjoint_not(noto, querypropmap(oors))  # [!#6,S]
    @test noto_.key === :or
    @test length(noto_.value) == 2
    @test noto_.value[1].key == :not
    @test noto_.value[2] == QueryLiteral(:symbol, :S)
    notarom = QueryOperator(:not, [QueryLiteral(:aromatic)])  # [A]
    notarom_ = resolve_disjoint_not(notarom, querypropmap(oors))  # [A]
    @test notarom_.key === :not
    @test notarom_.value[1] == QueryLiteral(:aromatic)
end

@testset "resolve_recursive" begin
    rec1 = QueryLiteral(:recursive, "[#6][NH]")  # [$([#6][NH])]
    rec2 = QueryLiteral(:recursive, "[#6]N")  # [$([#6]N)]
    rec1_ = resolve_recursive(rec1, querypropmap(rec2))  # [$([#6][NH]);#6;$([#6]N)]
    @test rec1_.key === :and
    @test rec1_.value[1] == QueryLiteral(:recursive, "[#6][NH]")
    @test rec1_.value[2] == QueryLiteral(:symbol, :C)
    @test rec1_.value[3] == QueryLiteral(:recursive, "[#6]N")
    rec2_ = resolve_recursive(rec2, querypropmap(rec1))  # [$([#6]N);#6]
    @test rec2_.key === :and
    @test rec2_.value[1] == QueryLiteral(:recursive, "[#6]N")
    @test rec2_.value[2] == QueryLiteral(:symbol, :C)
    aromn = QueryOperator(:and, [
        QueryLiteral(:symbol, :N),
        QueryLiteral(:isaromatic),
        QueryLiteral(:total_hydrogens, 1)
    ])  # [nH]
    nonan = QueryOperator(:and, [
        QueryLiteral(:symbol, :N),
        QueryOperator(:not, [QueryLiteral(:isaromatic)])
    ])  # N
    aromn_ = resolve_recursive(aromn, querypropmap(nonan))  # [nH]
    nonan_ = resolve_recursive(nonan, querypropmap(aromn))  # N
    @test aromn_.key === :and
    @test nonan_.value[2].key === :not
end

@testset "generate_truthtable" begin
    anytrue = QueryTree(QueryAny(true))
    tt1, tt2 = generate_truthtable(anytrue, anytrue)
    @test isempty(tt1.props)
    @test tt2.func([])
end

@testset "querymatch" begin
    op1 = QueryOperator(:and, [
        QueryLiteral(:symbol, :C),
        QueryOperator(:not, [QueryLiteral(:isaromatic)]),
        QueryOperator(:or, [
            QueryLiteral(:connectivity, 3),
            QueryLiteral(:connectivity, 4)
        ])
    ])
    @test length(querypropmap(op1)) == 3
    tbl1 = QueryTruthTable(op1)
    func1 = QueryTruthTable(
        v -> v[4] & ~v[3] & (v[1] | v[2]),
        [(:connectivity, 3), (:connectivity, 4), (:isaromatic,), (:symbol, :C)]
    )
    @test querymatch(tbl1, func1, true)
end

end