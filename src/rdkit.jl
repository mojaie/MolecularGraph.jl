#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    to_rdkmol, moltosmiles,
    morgan_fp_vector, rdkit_fp_vector,
    pattern_fp_vector, atom_pair_fp_vector,
    topological_torsion_fp_vector

using RDKitMinimalLib: Mol, get_mol, get_smiles,
    get_morgan_fp_as_bytes, get_rdkit_fp_as_bytes,
    get_pattern_fp_as_bytes, get_atom_pair_fp_as_bytes,
    get_topological_torsion_fp_as_bytes


function to_rdkmol(mol::MolGraph)
    return get_mol(to_json(Val{:rdkit}(), mol))
end


function moltosmiles(mol::MolGraph, details=nothing)
    return get_smiles(to_rdkmol(mol), details)
end


function uint8vec_to_bitarray(uvec::Vector{UInt8})
    bits = BitVector(undef, 8 * length(uvec))
    for (i, byte) in enumerate(uvec)
        for j in 0:7
            bits[8*(i-1) + j + 1] = (byte >> (7 - j)) & 0x01 == 1
        end
    end
    return bits
end


morgan_fp_vector(mol::Mol; kwargs...
    ) = uint8vec_to_bitarray(get_morgan_fp_as_bytes(mol::Mol, kwargs...))
morgan_fp_vector(mol::MolGraph; kwargs...
    ) = morgan_fp_vector(to_rdkmol(mol); kwargs...)

rdkit_fp_vector(mol::Mol; kwargs...
    ) = uint8vec_to_bitarray(get_rdkit_fp_as_bytes(mol::Mol, kwargs...))
rdkit_fp_vector(mol::MolGraph; kwargs...
    ) = rdkit_fp_vector(to_rdkmol(mol); kwargs...)

pattern_fp_vector(mol::Mol; kwargs...
    ) = uint8vec_to_bitarray(get_pattern_fp_as_bytes(mol::Mol, kwargs...))
pattern_fp_vector(mol::MolGraph; kwargs...
    ) = pattern_fp_vector(to_rdkmol(mol); kwargs...)

atom_pair_fp_vector(mol::Mol; kwargs...
    ) = uint8vec_to_bitarray(get_atom_pair_fp_as_bytes(mol::Mol, kwargs...))
atom_pair_fp_vector(mol::MolGraph; kwargs...
    ) = atom_pair_fp_vector(to_rdkmol(mol); kwargs...)

topological_torsion_fp_vector(mol::Mol; kwargs...
    ) = uint8vec_to_bitarray(get_topological_torsion_fp_as_bytes(mol::Mol, kwargs...))
topological_torsion_fp_vector(mol::MolGraph; kwargs...
    ) = topological_torsion_fp_vector(to_rdkmol(mol); kwargs...)


jaccard_index(a::BitVector, b::BitVector
    ) = sum(a .& b) / sum(a.| b)  # Not exported