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
    morgan_fp_vector(mol::RDKitMinimalLib.Mol; kwargs...) -> BitArray
    morgan_fp_vector(mol::SimpleMolGraph; kwargs...) -> BitArray

Return a Morgan fingerprint bit array
"""
MolecularGraph.morgan_fp_vector(mol::Mol; kwargs...
    ) = uint8vec_to_bitarray(get_morgan_fp_as_bytes(mol::Mol, kwargs...))
MolecularGraph.morgan_fp_vector(mol::SimpleMolGraph; kwargs...
    ) = morgan_fp_vector(rdkitmol(mol); kwargs...)


"""
    rdkit_fp_vector(mol::RDKitMinimalLib.Mol; kwargs...) -> BitArray
    rdkit_fp_vector(mol::SimpleMolGraph; kwargs...) -> BitArray

Return a RDKit fingerprint bit array
"""
MolecularGraph.rdkit_fp_vector(mol::Mol; kwargs...
    ) = uint8vec_to_bitarray(get_rdkit_fp_as_bytes(mol::Mol, kwargs...))
MolecularGraph.rdkit_fp_vector(mol::SimpleMolGraph; kwargs...
    ) = rdkit_fp_vector(rdkitmol(mol); kwargs...)


"""
    pattern_fp_vector(mol::RDKitMinimalLib.Mol; kwargs...) -> BitArray
    pattern_fp_vector(mol::SimpleMolGraph; kwargs...) -> BitArray

Return a pattern fingerprint bit array, a topological fingerprint
optimized for substructure screening
"""
MolecularGraph.pattern_fp_vector(mol::Mol; kwargs...
    ) = uint8vec_to_bitarray(get_pattern_fp_as_bytes(mol::Mol, kwargs...))
MolecularGraph.pattern_fp_vector(mol::SimpleMolGraph; kwargs...
    ) = pattern_fp_vector(rdkitmol(mol); kwargs...)


"""
    atom_pair_fp_vector(mol::RDKitMinimalLib.Mol; kwargs...) -> BitArray
    atom_pair_fp_vector(mol::SimpleMolGraph; kwargs...) -> BitArray

Return a atom pairs fingerprint bit array
"""
MolecularGraph.atom_pair_fp_vector(mol::Mol; kwargs...
    ) = uint8vec_to_bitarray(get_atom_pair_fp_as_bytes(mol::Mol, kwargs...))
MolecularGraph.atom_pair_fp_vector(mol::SimpleMolGraph; kwargs...
    ) = atom_pair_fp_vector(rdkitmol(mol); kwargs...)


"""
    topological_torsion_fp_vector(mol::RDKitMinimalLib.Mol; kwargs...) -> BitArray
    topological_torsion_fp_vector(mol::SimpleMolGraph; kwargs...) -> BitArray

Return a topological torsions fingerprint bit array
"""
MolecularGraph.topological_torsion_fp_vector(mol::Mol; kwargs...
    ) = uint8vec_to_bitarray(get_topological_torsion_fp_as_bytes(mol::Mol, kwargs...))
MolecularGraph.topological_torsion_fp_vector(mol::SimpleMolGraph; kwargs...
    ) = topological_torsion_fp_vector(rdkitmol(mol); kwargs...)


end # module