
@testset "svg" begin
    @testset "draw" begin
        ASSETS_DIR = joinpath(dirname(@__FILE__), "..", "..", "assets")
        demomol = sdftomol(open(joinpath(ASSETS_DIR, "test", "demo.mol")))
        drawsvg(demomol)
        # dest = open(joinpath(ASSETS_DIR, "image", "demo.svg"), "w")
        # write(dest, drawsvg(demomol))
    end
end
