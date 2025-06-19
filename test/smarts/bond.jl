
@testset "smarts.bond" begin

@testset "bondsymbol" begin
    state = SMARTSParser{Int,QueryAtom,QueryBond}("")
    qtree = QueryBond()
    implicit1 = bondsymbol!(state, qtree)
    @test implicit1 == 0

    state = SMARTSParser{Int,QueryAtom,QueryBond}("-")
    qtree = QueryBond()
    explicit1 = bondsymbol!(state, qtree)
    @test explicit1 == 1
    @test qtree == QueryBond(
        [(1, 2), (1, 3), (3, 4)],
        [qand(), qeq(:order, "1"), qnot(), qtrue(:isaromatic)])

    state = SMARTSParser{Int,QueryAtom,QueryBond}(raw"\?")
    qtree = QueryBond()
    stereo = bondsymbol!(state, qtree)
    @test explicit1 == 1
    @test qtree == QueryBond([(1, 2)], [qnot(), qeq(:stereo, "up")])
    @test state.pos == 3
end

@testset "smiles" begin
    state = SMILESParser{Int,SMILESAtom,SMILESBond}("#")
    triple = bond!(state)
    @test triple == SMILESBond(;order=3)

    state = SMILESParser{Int,SMILESAtom,SMILESBond}(":")
    arom = bond!(state)
    @test arom == SMILESBond(;isaromatic=true)

    state = SMILESParser{Int,SMILESAtom,SMILESBond}("\\")
    down = bond!(state)
    @test down == SMILESBond(;direction=:down)
end

@testset "smarts" begin
    state = SMARTSParser{Int,QueryAtom,QueryBond}("~")
    anyb = bond!(state)
    @test anyb == QueryBond(Tuple{Int,Int}[], [qanytrue()])

    state = SMARTSParser{Int,QueryAtom,QueryBond}("-!@")
    notring = bond!(state)
    @test notring == QueryBond(
        [(1, 2), (1, 3), (2, 4), (2, 5), (5, 6), (3, 7)],
        [qand(), qand(), qnot(), qeq(:order, "1"), qnot(),
        qtrue(:isaromatic), qtrue(:is_in_ring)])
end

end # smiles.bond
