
@testset "basedescriptor" begin

@testset "valence" begin
    alkyne = smilestomol("CC=CC#C")
    @test [a.pi for a in atomvector(alkyne)] == [0, 1, 1, 2, 2]
    pyrrole = smilestomol("C1=CC=CN1")
    @test [a.pi for a in atomvector(pyrrole)] == [1, 1, 1, 1, 0]
    azide = smilestomol("CC(=O)N=N=N")
    @test [a.pi for a in atomvector(azide)] == [0, 1, 1, 1, 2, 1]
    amide = smilestomol("CCC(=O)N")
    @test getatom(amide, 4).Hacceptor
    @test getatom(amide, 5).Hacceptor
    @test getatom(amide, 5).Hdonor
    fluoro = smilestomol("CCN(CO)CF")
    @test getatom(fluoro, 5).Hacceptor
    @test getatom(fluoro, 7).Hacceptor
    @test !getatom(fluoro, 3).Hdonor
    @test getatom(fluoro, 5).Hdonor
end

@testset "rotatable" begin
    Phe = smilestomol("N[C@@H](CC1=CC=CC=C1)C(O)=O")
    @test rotatable_count(Phe) == 3
    KCl = smilestomol("[K+].[Cl-]")
    @test rotatable_count(KCl) == 0
    dipyridamole = smilestomol(
        "n3c(nc2c(nc(nc2N1CCCCC1)N(CCO)CCO)c3N4CCCCC4)N(CCO)CCO")
    @test rotatable_count(dipyridamole) == 12
    paclitaxel = smilestomol(
        "CC1=C2[C@@]([C@]([C@H]([C@@H]3[C@]4([C@H](OC4)C[C@@H]([C@]3(C(=O)[C@@H]2OC(=O)C)C)O)OC(=O)C)OC(=O)c5ccccc5)(C[C@@H]1OC(=O)[C@H](O)[C@@H](NC(=O)c6ccccc6)c7ccccc7)O)(C)C"
    )
    @test rotatable_count(paclitaxel) == 15
end

end # basedescriptor
