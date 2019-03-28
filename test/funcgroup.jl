
using MolecularGraph: fgroupquery, fgroupcond, fgrouprecord

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
    fg = FunctionalGroup()
    gr = fg.nodeset
    gr[:Carbonyl] = fgroupquery(coumarin, "[#6]=[OD1]")
    res1 = fgrouprecord(
        coumarin, fg, Dict("have" => ["Carbonyl"], "query" => "[#6](=O)[#8]"))
    @test issetequal(collect(res1)[1], [3, 4, 5])

    gr[:Something] = Set{Set{Int}}()
    res2 = fgrouprecord(
        coumarin, fg,
        Dict("have" => ["Carbonyl", "Something"], "query" => "[#6](=O)[#8]"))
    @test isempty(res2)

    gr[:Carboxyl] = res1
    res3 = fgrouprecord(
        coumarin, fg, Dict("isa" => ["Carboxyl"], "query" => "c(=O)o"))
    @test issetequal(collect(res3)[1], [3, 4, 5])
    res4 = fgrouprecord(
        coumarin, fg, Dict("isa" => ["Carboxyl"], "query" => "[#6](=O)O"))
    @test isempty(res4)
end

@testset "functionalgroup" begin
    paclitaxel = smilestomol(
        "CC1=C2[C@@]([C@]([C@H]([C@@H]3[C@]4([C@H](OC4)C[C@@H]([C@]3(C(=O)[C@@H]2OC(=O)C)C)O)OC(=O)C)OC(=O)c5ccccc5)(C[C@@H]1OC(=O)[C@H](O)[C@@H](NC(=O)c6ccccc6)c7ccccc7)O)(C)C"
    )
    fg = functionalgroup(paclitaxel)
    res = fg.nodeset
    @test length(res[:Carboxyl]) == 4
    @test length(res[:SecAlcohol]) == 2
    @test length(res[:FourMemberedRing]) == 1
    @test length(res[:EightMemberedRing]) == 1
end

end # funcgroup
