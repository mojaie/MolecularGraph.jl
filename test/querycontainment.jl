
using MolecularGraph: resolve_disjoint_not, resolve_recursive, querypropmap

@testset "querycontainment" begin

@testset "resolve_disjoint_not" begin
    oors = QueryOperator(:or, [
        QueryLiteral(:symbol, :O), QueryLiteral(:symbol, :S)
    ])
    notn = QueryOperator(:not, [QueryLiteral(:symbol, :N)])
    oors_ = resolve_disjoint_not(oors, querypropmap(notn))
    notn_ = resolve_disjoint_not(notn, querypropmap(oors))
    @test notn_.key === :or
    @test notn_.value[1].key == :not
    @test notn_.value[1].value[1] == QueryLiteral(:symbol, :N)
    @test notn_.value[2] == QueryLiteral(:symbol, :O)
    @test notn_.value[3] == QueryLiteral(:symbol, :S)
    @test oors_.key === :or
    @test oors_.value[1] == QueryLiteral(:symbol, :O)
    @test oors_.value[2] == QueryLiteral(:symbol, :S)
    noto = QueryOperator(:not, [QueryLiteral(:symbol, :O)])
    noto_ = resolve_disjoint_not(noto, querypropmap(oors))
    @test noto_.key === :or
    @test length(noto_.value) == 2
    @test noto_.value[1].key == :not
    @test noto_.value[2] == QueryLiteral(:symbol, :S)
    notarom = QueryOperator(:not, [QueryLiteral(:aromatic)])
    notarom_ = resolve_disjoint_not(notarom, querypropmap(oors))
    @test notarom_.key === :not
    @test notarom_.value[1] == QueryLiteral(:aromatic)
end

@testset "resolve_recursive" begin
    rec1 = QueryLiteral(:recursive, "[#6][NH]")
    rec2 = QueryLiteral(:recursive, "[#6]N")
    rec1_ = resolve_recursive(rec1, querypropmap(rec2))
    rec2_ = resolve_recursive(rec2, querypropmap(rec1))
    @test rec1_.key === :and
    @test rec1_.value[1] == QueryLiteral(:recursive, "[#6][NH]")
    @test rec1_.value[2] == QueryLiteral(:symbol, :C)
    @test rec1_.value[3] == QueryLiteral(:recursive, "[#6]N")
    @test rec2_.key === :and
    @test rec2_.value[1] == QueryLiteral(:recursive, "[#6]N")
    @test rec2_.value[2] == QueryLiteral(:symbol, :C)
    aromn = QueryOperator(:and, [
        QueryLiteral(:symbol, :N),
        QueryLiteral(:isaromatic),
        QueryLiteral(:total_hydrogens, 1)
    ])
    nonan = QueryOperator(:and, [
        QueryLiteral(:symbol, :N),
        QueryOperator(:not, [QueryLiteral(:isaromatic)])
    ])
    aromn_ = resolve_recursive(aromn, querypropmap(nonan))
    nonan_ = resolve_recursive(nonan, querypropmap(aromn))
    @test aromn_.key === :and
    @test nonan_.value[2].key === :not
end

end