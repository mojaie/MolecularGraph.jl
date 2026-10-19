
@testset "draw3d" begin
    mol = sdftomol(joinpath(@__DIR__, "..", "..", "assets", "test", "aspirin_3d.sdf"))
    @test length(atom_radius(mol)) == nv(mol)
    @test atom_radius(mol; mapping=0.5) == fill(0.5, nv(mol))
    # carbon radius by hybridization, smallest value for metals with multiple values
    m = smilestomol("CC=CC#C[Fe]")
    @test atom_radius(m; mapping=ATOM_COVALENT_RADII) ≈ [0.76, 0.73, 0.73, 0.69, 0.69, 1.32] atol=1e-6
    @test atom_radius(smilestomol("CO")) ≈ [ATOM_VANDERWAALS_RADII[6], ATOM_VANDERWAALS_RADII[8]]

    # MakieExt (requires Makie, no backend is needed to build plots)
    for f in (spacefilling, ballstick, stick, wire)
        @test f(mol) isa Makie.FigureAxisPlot
    end
    fig = Makie.Figure()
    ls = Makie.LScene(fig[1, 1])
    @test length(spacefilling!(ls, mol).plots) == 1  # atoms only
    for f in (ballstick!, stick!, wire!)
        @test length(f(ls, mol).plots) > 1  # atoms and bonds
    end
    @test ballstick!(ls, mol; radii=0.5).radii[] == 0.5
    @test wire!(ls, mol; size=0.3).bonddiameter[] == 0.3
end
