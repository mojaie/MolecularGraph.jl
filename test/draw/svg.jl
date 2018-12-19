
@testset "svg" begin
    @testset "draw" begin
        ASSETS_DIR = joinpath(dirname(@__FILE__), "..", "..", "assets")
        demomol = sdftomol(open(joinpath(ASSETS_DIR, "test", "demo.mol")))
        drawsvg!(demomol, 200, 200)
        # dest = open(joinpath(ASSETS_DIR, "image", "demo.svg"), "w")
        # write(dest, drawsvg!(demomol, 200, 200))
    end
end
