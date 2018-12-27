
@testset "smarts.atom" begin

@testset "atomsymbol" begin
    state = ConnectedSmarts("Cl")
    Chloride = atomsymbol!(state)
    @test Chloride == (:and => (:Symbol => :Cl, :Aromatic => false))
    @test state.pos == 3

    state = ConnectedSmarts("Cr")
    Chromium = atomsymbol!(state)
    @test Chromium == (:and => (:Symbol => :C, :Aromatic => false))
    @test state.pos == 2

    state = ConnectedSmarts("p")
    aromp = atomsymbol!(state)
    @test aromp == (:and => (:Symbol => :P, :Aromatic => true))
end

@testset "atomprop" begin
    state = ConnectedSmarts("s")
    aroms = atomprop!(state)
    @test aroms == (:and => (:Symbol => :S, :Aromatic => true))

    state = ConnectedSmarts("123")
    iso = atomprop!(state)
    @test iso == (:Mass => 123)

    state = ConnectedSmarts("2H]")
    forward!(state)
    hyd = atomprop!(state)
    @test hyd == (:Symbol => :H)

    state = ConnectedSmarts("!H]")
    forward!(state)
    noth = atomprop!(state)
    @test noth == (:H_Count => 1)

    state = ConnectedSmarts("H+")
    proton = atomprop!(state)
    @test proton == (:and => (:Symbol => :H, :Charge => 1))

    state = ConnectedSmarts("H41")
    H4 = atomprop!(state)
    @test H4 == (:H_Count => 4)
    @test state.pos == 3

    state = ConnectedSmarts("X2")
    X2 = atomprop!(state)
    @test X2 == (:Connectivity => 2)

    state = ConnectedSmarts("Xe")
    Xe = atomprop!(state)
    @test Xe == (:Symbol => :Xe)
    @test state.pos == 3

    state = ConnectedSmarts("Na")
    Na = atomprop!(state)
    @test Na == (:Symbol => :Na)
    @test state.pos == 3

    state = ConnectedSmarts("Yv2")
    Yval = atomprop!(state)
    @test Yval == (:Symbol => :Y)

    state = ConnectedSmarts("+23")
    chg1 = atomprop!(state)
    @test chg1 == (:Charge => 2)
    @test state.pos == 3

    state = ConnectedSmarts("----+")
    chg2 = atomprop!(state)
    @test chg2 == (:Charge => -4)
    @test state.pos == 5

    state = ConnectedSmarts("@+")
    stereo1 = atomprop!(state)
    @test stereo1 == (:stereo => 1)
    @test state.pos == 2

    state = ConnectedSmarts("@@?")
    stereo4 = atomprop!(state)
    @test stereo4 == (:stereo => 4)
    @test state.pos == 4

    state = ConnectedSmarts("\$([CH2]=*)")
    rec = atomprop!(state)
    @test rec == (:recursive => "[CH2]=*")
    @test state.pos == 11
end

@testset "atom" begin
    state = SmilesParser("Br")
    br = atom!(state)
    @test br.symbol == :Br

    state = SmilesParser("[14c@@H]")
    iso = atom!(state)
    @test iso.symbol == :C
    @test iso.mass == 14
    @test iso.isaromatic
    @test iso.stereo == 2
    @test state.pos == 9

    state = SmilesParser("[Zn++]")
    zn = atom!(state)
    @test zn.symbol == :Zn
    @test zn.charge == 2

    state = SmilesParser("[O-2]")
    ox = atom!(state)
    @test ox.symbol == :O
    @test ox.charge == -2
end

@testset "smartsatom" begin
    state = ConnectedSmarts("c")
    aromc = atom!(state)
    @test aromc.query == (:and => (:Symbol => :C, :Aromatic => true))

    state = ConnectedSmarts("[#16]")
    no16 = atom!(state)
    @test no16.query == (:Symbol => :S)
end

end # smiles.atom
