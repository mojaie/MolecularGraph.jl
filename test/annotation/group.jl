
@testset "annotation.group" begin

@testset "group" begin
    assign_annot!(mol) = (atom_annot!(mol); group_annot!(mol))
    nullmol = smilestomol("")
    assign_annot!(nullmol)
    @test count(nullmol.v[:Carbonyl]) == 0
    alcohol = smilestomol("CC(O)C(CO)(C)O")
    assign_annot!(alcohol)
    @test alcohol.v[:Oxygen] == [
        nothing, nothing, 1, nothing, nothing, 1, nothing, 1]
    @test alcohol.v[:Alcohol] == [
        nothing, 2, nothing, 3, 1, nothing, nothing, nothing]
    amine = smilestomol("CNC(C#N)(C=N)N")
    assign_annot!(amine)
    @test amine.v[:Nitrogen] == [
        nothing, 2, nothing, nothing, 8, nothing, 5, 1]
end

end # annotation.group
