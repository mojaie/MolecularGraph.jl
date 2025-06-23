#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

module RDKitExt

using MolecularGraph:
    MolecularGraph, SimpleMolGraph,
    to_json, rdkitmol, smiles,
    morgan_fp_vector, rdkit_fp_vector,
    pattern_fp_vector, atom_pair_fp_vector,
    topological_torsion_fp_vector

using RDKitMinimalLib:
    RDKitMinimalLib, Mol, get_mol, get_smiles,
    get_morgan_fp_as_bytes, get_rdkit_fp_as_bytes,
    get_pattern_fp_as_bytes, get_atom_pair_fp_as_bytes,
    get_topological_torsion_fp_as_bytes


function uint8vec_to_bitarray(uvec::Vector{UInt8})
    bits = BitVector(undef, 8 * length(uvec))
    for (i, byte) in enumerate(uvec)
        for j in 1:8
            bits[8*(i-1) + j] = (byte >> (8 - j)) & 0x01
        end
    end
    return bits
end


"""
    rdkitmol(mol::SimpleMolGraph) -> RDKitMinimalLib.Mol

Convert the molecule object into a RDKit molecule object that can be used in RDKitMinimalLib.jl
"""
function MolecularGraph.rdkitmol(mol::SimpleMolGraph)
    return get_mol(to_json(Val{:rdkit}(), mol))
end


"""
    smiles(mol::SimpleMolGraph) -> String

Return a SMILES string of the molecule object.
"""
function MolecularGraph.smiles(mol::SimpleMolGraph, details=nothing)
    return get_smiles(rdkitmol(mol), details)
end


"""
    morgan_fp_vector(mol::RDKitMinimalLib.Mol, details=nothing) -> BitArray
    morgan_fp_vector(mol::SimpleMolGraph, details=nothing) -> BitArray

Return a Morgan fingerprint bit array
"""
MolecularGraph.morgan_fp_vector(mol::Mol, details=nothing
    ) = uint8vec_to_bitarray(get_morgan_fp_as_bytes(mol::Mol, details))
MolecularGraph.morgan_fp_vector(mol::SimpleMolGraph, details=nothing
    ) = morgan_fp_vector(rdkitmol(mol), details)


"""
    rdkit_fp_vector(mol::RDKitMinimalLib.Mol, details=nothing) -> BitArray
    rdkit_fp_vector(mol::SimpleMolGraph, details=nothing) -> BitArray

Return a RDKit fingerprint bit array
"""
MolecularGraph.rdkit_fp_vector(mol::Mol, details=nothing
    ) = uint8vec_to_bitarray(get_rdkit_fp_as_bytes(mol::Mol, details))
MolecularGraph.rdkit_fp_vector(mol::SimpleMolGraph, details=nothing
    ) = rdkit_fp_vector(rdkitmol(mol), details)


"""
    pattern_fp_vector(mol::RDKitMinimalLib.Mol, details=nothing) -> BitArray
    pattern_fp_vector(mol::SimpleMolGraph, details=nothing) -> BitArray

Return a pattern fingerprint bit array, a topological fingerprint
optimized for substructure screening
"""
MolecularGraph.pattern_fp_vector(mol::Mol, details=nothing
    ) = uint8vec_to_bitarray(get_pattern_fp_as_bytes(mol::Mol, details))
MolecularGraph.pattern_fp_vector(mol::SimpleMolGraph, details=nothing
    ) = pattern_fp_vector(rdkitmol(mol), details)


"""
    atom_pair_fp_vector(mol::RDKitMinimalLib.Mol, details=nothing) -> BitArray
    atom_pair_fp_vector(mol::SimpleMolGraph, details=nothing) -> BitArray

Return a atom pairs fingerprint bit array
"""
MolecularGraph.atom_pair_fp_vector(mol::Mol, details=nothing
    ) = uint8vec_to_bitarray(get_atom_pair_fp_as_bytes(mol::Mol, details))
MolecularGraph.atom_pair_fp_vector(mol::SimpleMolGraph, details=nothing
    ) = atom_pair_fp_vector(rdkitmol(mol), details)


"""
    topological_torsion_fp_vector(mol::RDKitMinimalLib.Mol, details=nothing) -> BitArray
    topological_torsion_fp_vector(mol::SimpleMolGraph, details=nothing) -> BitArray

Return a topological torsions fingerprint bit array
"""
MolecularGraph.topological_torsion_fp_vector(mol::Mol, details=nothing
    ) = uint8vec_to_bitarray(get_topological_torsion_fp_as_bytes(mol::Mol, details))
MolecularGraph.topological_torsion_fp_vector(mol::SimpleMolGraph, details=nothing
    ) = topological_torsion_fp_vector(rdkitmol(mol), details)


end # module