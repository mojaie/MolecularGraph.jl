
@testset "svg" begin
    @testset "draw" begin
        ASSET_DIR = joinpath(dirname(@__FILE__), "..", "..", "assets", "test")
        demomol = loadsdfmol(open(joinpath(ASSET_DIR, "demo.mol")))
        print(drawsvg!(demomol, 200, 200))
    end
end
