
@testset "svg" begin
    @testset "draw" begin
        ASSETS_DIR = joinpath(dirname(@__FILE__), "..", "..", "assets")
        demomol = sdftomol(open(joinpath(ASSETS_DIR, "test", "demo.mol")))
        drawsvg(demomol)
        nullmol = sdftomol(open(joinpath(ASSETS_DIR, "test", "null.mol")))
        drawsvg(nullmol)
        # dest = open(joinpath(ASSETS_DIR, "image", "demo.svg"), "w")
        # write(dest, drawsvg(demomol))
    end
end
