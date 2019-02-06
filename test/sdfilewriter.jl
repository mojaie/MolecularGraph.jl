
@testset "sdfilewriter" begin

@testset "molblock" begin
    demomol = joinpath(dirname(@__FILE__), "..", "assets", "test", "demo.mol")
    mol = sdftomol(demomol)
    buf = IOBuffer(write=true)
    sdfilewriter(buf, (mol,))
    mol2 = sdftomol(split(String(take!(buf)), "\n"))
    @test nodecount(mol2) == 37
end

end #sdfilewriter
