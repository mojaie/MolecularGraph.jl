#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

# RDKitExt (requires RDKitMinimalLib), migrated from ext/tests_ext.jl

@testset "ext.rdkit" begin
    ASSETS_DIR = joinpath(dirname(@__FILE__), "..", "..", "assets")
    mol1 = sdftomol(open(joinpath(ASSETS_DIR, "test", "nata.mol")))
    # Tacrolimus
    mol2 = smilestomol(raw"C=CC[C@@H]1C=C(C)C[C@H](C)C[C@H](OC)[C@H]2O[C@](O)([C@H](C)C[C@@H]2OC)C(=O)C(=O)N2CCCC[C@H]2C(=O)O[C@H](\C(=C\[C@@H]2CC[C@@H](O)[C@H](OC)C2)/C)[C@H](C)[C@@H](O)CC1=O")

    @testset "conversion" begin
        @test rdkitmol(mol1) isa RDKitMinimalLib.Mol
        @test is_ring_aromatic(mol_from_json(JSON.json(to_rdkdict(mol1)))) == is_ring_aromatic(mol1)
        @test is_ring_aromatic(smilestomol(smiles(mol1))) == is_ring_aromatic(mol1)
        aspirin = smilestomol("CC(=O)Oc1ccccc1C(=O)O")
        @test smiles(aspirin) == "CC(=O)Oc1ccccc1C(=O)O"
    end

    @testset "fingerprint" begin
        @test length(morgan_fp_string(mol1)) == 2048
        @test length(morgan_fp_string(mol2)) == 2048
        a = morgan_fp_vector(mol1)
        b = morgan_fp_vector(mol2)
        @test sum(a .& b) == 21
        @test sum(a .| b) == 144
        a = rdkit_fp_vector(mol1)
        b = rdkit_fp_vector(mol2)
        @test sum(a .& b) == 1497
        @test sum(a .| b) == 1920
        a = pattern_fp_vector(mol1)
        b = pattern_fp_vector(mol2)
        @test sum(a .& b) == 1058
        @test sum(a .| b) == 1240
        # *_fp_vector (from bytes, bit-reversed) should be consistent with *_fp_string
        for (vecf, strf) in [
                (morgan_fp_vector, morgan_fp_string),
                (rdkit_fp_vector, rdkit_fp_string),
                (pattern_fp_vector, pattern_fp_string),
                (atom_pair_fp_vector, atom_pair_fp_string),
                (topological_torsion_fp_vector, topological_torsion_fp_string)]
            @test strf(mol2) == join(Int.(vecf(mol2)))
        end
    end
end
