
@testset "annotation.elemental" begin

@testset "elemental" begin
    alkyne = smilestomol("CC=CC#C")
    @test alkyne.v[:Pi] == [0, 1, 1, 2, 2]
    pyrrole = smilestomol("C1=CC=CN1")
    @test pyrrole.v[:Pi] == [1, 1, 1, 1, 0]
    azide = smilestomol("CC(=O)N=N=N")
    @test azide.v[:Pi] == [0, 1, 1, 1, 2, 1]
    amide = smilestomol("CCC(=O)N")
    @test amide.v[:H_Acceptor] == [0, 0, 0, 1, 1]
    @test amide.v[:H_Donor] == [0, 0, 0, 0, 1]
    fluoro = smilestomol("CCN(CO)CF")
    @test fluoro.v[:H_Acceptor] == [0, 0, 1, 0, 1, 0, 1]
    @test fluoro.v[:H_Donor] == [0, 0, 0, 0, 1, 0, 0]
end

end # annotation.elemental
