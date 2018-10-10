
@testset "svg" begin
    ASSET_DIR = joinpath(dirname(@__FILE__), "..", "..", "assets", "test")
    demomol = loadsdfmol(open(joinpath(ASSET_DIR, "demo.mol")))
    draw(:svg, demomol)
end
