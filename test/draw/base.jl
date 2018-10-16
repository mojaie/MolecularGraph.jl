
@testset "draw.base" begin
    @testset "base" begin
        ASSET_DIR = joinpath(dirname(@__FILE__), "..", "..", "assets", "test")
        demomol = loadsdfmol(open(joinpath(ASSET_DIR, "demo.mol")))
        display_terminal_carbon!(demomol)
        atoms = atomvector(demomol)
        @test sum(1 for a in atoms if a.symbol == "C" && a.visible) == 3
        equalize_terminal_double_bond!(demomol)
        bonds = bondvector(demomol)
        @test sum(1 for b in bonds if b.order == 2 && b.notation == 2) == 2
        double_bond_along_ring!(demomol)
        @test sum(1 for b in bonds if b.order == 2 && b.notation == 1) == 1
    end
end
