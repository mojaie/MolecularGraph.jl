
@testset "svg" begin
    @testset "draw" begi
        ASSET_DIR = joinpath(dirname(@__FILE__), "..", "..", "assets", "test")
        demomol = loadsdfmol(open(joinpath(ASSET_DIR, "demo.mol")))
        draw(:svg, demomol)
    end
end
