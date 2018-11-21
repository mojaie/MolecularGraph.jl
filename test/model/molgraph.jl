
@testset "molgraph" begin

    @testset "mapmol" begin
        sdf = GMapMol{SDFileAtom,SDFileBond}()
        vmol = vectormol(sdf)
        @test vmol isa GVectorMol{SDFileAtom,SDFileBond}
    end
end
