
@testset "svg" begin
    @testset "draw" begin
        IMG_DIR = joinpath(dirname(@__FILE__), "..", "..", "assets", "image")
        demomol = loadsdfmol(open(joinpath(IMG_DIR, "demo.mol")))
        # dest = open(joinpath(IMG_DIR, "demo.svg"), "w")
        # write(dest, drawsvg!(demomol, 200, 200))
    end
end
