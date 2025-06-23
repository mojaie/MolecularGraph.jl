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
    get_morgan_fp, get_morgan_fp_as_bytes,
    get_rdkit_fp, get_rdkit_fp_as_bytes,
    get_pattern_fp, get_pattern_fp_as_bytes,
    get_atom_pair_fp, get_atom_pair_fp_as_bytes,
    get_topological_torsion_fp, get_topological_torsion_fp_as_bytes


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
    RDKitMinimalLib.get_smiles(mol::SimpleMolGraph) -> String
    smiles(mol::SimpleMolGraph) -> String

Return a SMILES string of the molecule object.
"""
RDKitMinimalLib.get_smiles(mol::SimpleMolGraph, args...
    ) = get_smiles(rdkitmol(mol), args...)
MolecularGraph.smiles(mol, args...) = RDKitMinimalLib.get_smiles(mol, args...)


"""
    morgan_fp_vector(mol::RDKitMinimalLib.Mol, details=nothing) -> BitArray
    morgan_fp_vector(mol::SimpleMolGraph, details=nothing) -> BitArray

Return a Morgan fingerprint bit array
"""
MolecularGraph.morgan_fp_vector(mol::Mol, args...
    ) = uint8vec_to_bitarray(get_morgan_fp_as_bytes(mol::Mol, args...))
MolecularGraph.morgan_fp_vector(mol::SimpleMolGraph, args...
    ) = morgan_fp_vector(rdkitmol(mol), args...)

"""
    RDKitMinimalLib.get_morgan_fp(mol::SimpleMolGraph, details=nothing) -> String
    morgan_fp_string(mol::SimpleMolGraph, details=nothing) -> String

Return a Morgan fingerprint string
"""
RDKitMinimalLib.get_morgan_fp(mol::SimpleMolGraph, args...
    ) = get_morgan_fp(rdkitmol(mol), args...)
MolecularGraph.morgan_fp_string(args...) = get_morgan_fp(args...)


"""
    rdkit_fp_vector(mol::RDKitMinimalLib.Mol, details=nothing) -> BitArray
    rdkit_fp_vector(mol::SimpleMolGraph, details=nothing) -> BitArray

Return a RDKit fingerprint bit array
"""
MolecularGraph.rdkit_fp_vector(mol::Mol, args...
    ) = uint8vec_to_bitarray(get_rdkit_fp_as_bytes(mol::Mol, args...))
MolecularGraph.rdkit_fp_vector(mol::SimpleMolGraph, args...
    ) = rdkit_fp_vector(rdkitmol(mol), args...)

"""
    RDKitMinimalLib.get_rdkit_fp(mol::SimpleMolGraph, details=nothing) -> String
    rdkit_fp_string(mol::SimpleMolGraph, details=nothing) -> String

Return a rdkit fingerprint string
"""
RDKitMinimalLib.get_rdkit_fp(mol::SimpleMolGraph, args...
    ) = get_rdkit_fp(rdkitmol(mol), args...)
MolecularGraph.rdkit_fp_string(args...) = get_rdkit_fp(args...)


"""
    pattern_fp_vector(mol::RDKitMinimalLib.Mol, details=nothing) -> BitArray
    pattern_fp_vector(mol::SimpleMolGraph, details=nothing) -> BitArray

Return a pattern fingerprint bit array
"""
MolecularGraph.pattern_fp_vector(mol::Mol, args...
    ) = uint8vec_to_bitarray(get_pattern_fp_as_bytes(mol::Mol, args...))
MolecularGraph.pattern_fp_vector(mol::SimpleMolGraph, args...
    ) = pattern_fp_vector(rdkitmol(mol), args...)

"""
    RDKitMinimalLib.get_pattern_fp(mol::SimpleMolGraph, details=nothing) -> String
    pattern_fp_string(mol::SimpleMolGraph, details=nothing) -> String

Return a pattern fingerprint string
"""
RDKitMinimalLib.get_pattern_fp(mol::SimpleMolGraph, args...
    ) = get_pattern_fp(rdkitmol(mol), args...)
MolecularGraph.pattern_fp_string(args...) = get_pattern_fp(args...)


"""
    atom_pair_fp_vector(mol::RDKitMinimalLib.Mol, details=nothing) -> BitArray
    atom_pair_fp_vector(mol::SimpleMolGraph, details=nothing) -> BitArray

Return a atom pairs fingerprint bit array
"""
MolecularGraph.atom_pair_fp_vector(mol::Mol, args...
    ) = uint8vec_to_bitarray(get_atom_pair_fp_as_bytes(mol::Mol, args...))
MolecularGraph.atom_pair_fp_vector(mol::SimpleMolGraph, args...
    ) = atom_pair_fp_vector(rdkitmol(mol), args...)


"""
    RDKitMinimalLib.get_atom_pair_fp(mol::SimpleMolGraph, details=nothing) -> String
    atom_pair_fp_string(mol::SimpleMolGraph, details=nothing) -> String

Return a atom pairs fingerprint string
"""
RDKitMinimalLib.get_atom_pair_fp(mol::SimpleMolGraph, args...
    ) = get_atom_pair_fp(rdkitmol(mol), args...)
MolecularGraph.atom_pair_fp_string(args...) = get_atom_pair_fp(args...)





"""
    topological_torsion_fp_vector(mol::RDKitMinimalLib.Mol, details=nothing) -> BitArray
    topological_torsion_fp_vector(mol::SimpleMolGraph, details=nothing) -> BitArray

Return a atom pairs fingerprint bit array
"""
MolecularGraph.topological_torsion_fp_vector(mol::Mol, args...
    ) = uint8vec_to_bitarray(get_topological_torsion_fp_as_bytes(mol::Mol, args...))
MolecularGraph.topological_torsion_fp_vector(mol::SimpleMolGraph, args...
    ) = topological_torsion_fp_vector(rdkitmol(mol), args...)


"""
    RDKitMinimalLib.get_topological_torsion_fp(mol::SimpleMolGraph, details=nothing) -> String
    topological_torsion_fp_string(mol::SimpleMolGraph, details=nothing) -> String

Return a atom pairs fingerprint string
"""
RDKitMinimalLib.get_topological_torsion_fp(mol::SimpleMolGraph, args...
    ) = get_topological_torsion_fp(rdkitmol(mol), args...)
MolecularGraph.topological_torsion_fp_string(args...) = get_topological_torsion_fp(args...)


end # module