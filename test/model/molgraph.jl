
@testset "molgraph" begin

    @testset "mapmol" begin
        sdf = mapmol(SDFileAtom, SDFileBond)
        vmol = vectormol(sdf)
        @test vmol isa VectorMol{SDFileAtom,SDFileBond}
    end

    @testset "vectormol" begin
        mol = vectormol(SDFileAtom, SDFileBond)
        @test nodecount(mol) == 0
        atoms = [SDFileAtom(),SDFileAtom(),SDFileAtom()]
        bonds = [SDFileBond(1, 2)]
        mol = vectormol(atoms, bonds)
        newbond = setnodes(getedge(mol, 1), 2, 3)
        # TODO: updateedge! for vectorgraph
        # updateedge!(mol, newbond, 2)
        # @test edgecount(mol) == 2
    end

    @testset "querymol" begin
        q = SMARTS()
        a = SmartsAtom(:Any => true)
        a2 = SmartsAtom(:Any => true)
        b = SmartsBond(1, 2, :Any => true)
        updateatom!(q, a, 1)
        updateatom!(q, a2, 2)
        updatebond!(q, b, 1)
        @test atomcount(q) == 2
    end

end
