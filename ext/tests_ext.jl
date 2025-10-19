#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

module MolecularGraphExtTest

using MolecularGraph
using Test

import JSON
import RDKitMinimalLib


@testset "rdkit" begin
    ASSETS_DIR = joinpath(dirname(@__FILE__), "..", "assets")
    mol1 = sdftomol(open(joinpath(ASSETS_DIR, "test", "nata.mol")))
    # Tacrolimus
    mol2 = smilestomol(raw"C=CC[C@@H]1C=C(C)C[C@H](C)C[C@H](OC)[C@H]2O[C@](O)([C@H](C)C[C@@H]2OC)C(=O)C(=O)N2CCCC[C@H]2C(=O)O[C@H](\C(=C\[C@@H]2CC[C@@H](O)[C@H](OC)C2)/C)[C@H](C)[C@@H](O)CC1=O")
    @test is_ring_aromatic(mol_from_json(JSON.json(to_rdkdict(mol1)))) == is_ring_aromatic(mol1)
    @test is_ring_aromatic(smilestomol(smiles(mol1))) == is_ring_aromatic(mol1)
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
end

end
