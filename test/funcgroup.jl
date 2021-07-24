
using MolecularGraph:
    fgroupquery, fgroupcond, fgrouprecord, getterm, setterm!

@testset "funcgroup" begin

@testset "fgroupquery" begin
    glyc = smilestomol("C(O)C(O)CO")
    res = fgroupquery(glyc, "[#6][OD1]")
    @test length(res) == 3
end

@testset "fgroupcond" begin
    alkoxide = smilestomol("C(O)C(O)C[O-]")
    res = fgroupcond(alkoxide, Dict("any" => ["[#6][OD1]", "[#6][O-]"]))
    @test length(res) == 3
end

@testset "fgrouprecord" begin
    coumarin = smilestomol("C1=CC(=O)OC2=C1C=CC=C2")
    fgc = FunctionalGroupClassifier()
    setterm!(fgc, :Carbonyl, fgroupquery(coumarin, "[#6]=[OD1]"))
    res1 = fgrouprecord(
        coumarin, fgc, Dict("has" => ["Carbonyl"], "query" => "[#6](=O)[#8]"))
    @test issetequal(collect(res1)[1], [3, 4, 5])

    res2 = fgrouprecord(
        coumarin, fgc,
        Dict("has" => ["Carbonyl", "Something"], "query" => "[#6](=O)[#8]"))
    @test isempty(res2)

    setterm!(fgc, :Carboxyl, res1)
    res3 = fgrouprecord(
        coumarin, fgc, Dict("isa" => ["Carboxyl"], "query" => "c(=O)o"))
    @test issetequal(collect(res3)[1], [3, 4, 5])
    res4 = fgrouprecord(
        coumarin, fgc, Dict("isa" => ["Carboxyl"], "query" => "[#6](=O)O"))
    @test isempty(res4)
end

end # funcgroup
