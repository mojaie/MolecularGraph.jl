
@testset "molgraph" begin

    @testset "graphmol" begin
        mol = graphmol(SDFileAtom, SDFileBond)
        @test nodecount(mol) == 0
        edges = [(1, 2)]
        atoms = [SDFileAtom(),SDFileAtom(),SDFileAtom()]
        bonds = [SDFileBond()]
        mol = graphmol(edges, atoms, bonds)
        newbond = addedge!(mol, 2, 3, SDFileBond())
        @test newbond == 2
        # TODO: updateedge! for plaingraph
        # updateedge!(mol, newbond, 2)
        # @test edgecount(mol) == 2
    end

    @testset "querymol" begin
        q = querymol(SmartsAtom, SmartsBond)
        addnode!(q, SmartsAtom(:Any => true))
        addnode!(q, SmartsAtom(:Any => true))
        addedge!(q, 1, 2, SmartsBond(:Any => true))
        @test nodecount(q) == 2
        @test hasedge(q, 1, 2)
    end

end
