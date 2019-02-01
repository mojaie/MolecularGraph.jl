
@testset "annotation.rotatable" begin

@testset "rotatable" begin
    Phe = smilestomol("N[C@@H](CC1=CC=CC=C1)C(O)=O")
    rotatable!(Phe)
    @test count(Phe[:Rotatable]) == 3
    KCl = smilestomol("[K+].[Cl-]")
    rotatable!(KCl)
    @test count(KCl[:Rotatable]) == 0
    dipyridamole = smilestomol(
        "n3c(nc2c(nc(nc2N1CCCCC1)N(CCO)CCO)c3N4CCCCC4)N(CCO)CCO")
    rotatable!(dipyridamole)
    @test count(dipyridamole[:Rotatable]) == 12
    paclitaxel = smilestomol(
        "CC1=C2[C@@]([C@]([C@H]([C@@H]3[C@]4([C@H](OC4)C[C@@H]([C@]3(C(=O)[C@@H]2OC(=O)C)C)O)OC(=O)C)OC(=O)c5ccccc5)(C[C@@H]1OC(=O)[C@H](O)[C@@H](NC(=O)c6ccccc6)c7ccccc7)O)(C)C"
    )
    rotatable!(paclitaxel)
    @test count(paclitaxel[:Rotatable]) == 15
end


end # annotation.rotatable
