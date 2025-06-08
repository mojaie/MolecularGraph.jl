#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

function rdk_on_init!(mol)
    set_state!(mol, :initialized, true)
end


function rdk_on_update!(mol)
    update_edge_rank!(mol)
    reset_updates!(mol)
    # Preprocess
    default_atom_charge!(mol)
    default_bond_order!(mol)
    # Cache relatively expensive descriptors
    sssr!(mol)
    apparent_valence!(mol)
    valence!(mol)
    lone_pair!(mol)
    is_ring_aromatic!(mol)
end


function molgraph_from_dict(
        ::Val{:rdkit}, ::Type{T}, ::Type{V}, ::Type{E}, data::Dict;
        on_init=rdk_on_init!, on_update=rdk_on_update!, kwargs...) where {T,V,E}
    gps = MolGraphProperty{T}()
    # edges
    es = Edge{T}[]
    eps = Dict{Edge{T},E}()
    iscis = Dict("cis" => true, "trans" => false)
    for b in data["bonds"]
        edge = u_edge(T, b["atoms"][1] + 1, b["atoms"][2] + 1)
        push!(es, edge)
        if haskey(b, "stereo") && b["stereo"] in keys(iscis)
            u, v = b["stereoAtoms"]
            gps.stereobond[edge] = (u + 1, v + 1, iscis[b["stereo"]])
        end
        eps[edge] = CommonChemBond(b)
    end
    g = SimpleGraph(es)
    # vertices
    vps = Dict{T,V}()
    isclockwise = Dict("cw" => true, "ccw" => false)
    for (i, a) in enumerate(data["atoms"])
        vps[i] = CommonChemAtom(a)
        if haskey(a, "stereo") && a["stereo"] in keys(isclockwise)
            nbrs = ordered_neighbors(g, i)
            # TODO: ambiguity in implicit H
            gps.stereocenter[i] = (nbrs[1:3]..., isclockwise[a["stereo"]])
        end
    end
    if haskey(data, "conformers")
        for cds in data["conformers"]
            # Wedges cannot be preserved. Use coordgen to generate 2D coords manually
            cds["dim"] == 3 || continue
            push!(gps.coords3d, [Point3d(cd...) for cd in cds["coords"]])
        end
    end
    return MolGraph{T,V,E}(
        g, vps, eps, gprops=gps, on_init=on_init, on_update=on_update; kwargs...)
end


function MolGraph{T,CommonChemAtom,CommonChemBond}(data::Dict; kwargs...) where T<:Integer
    if data["commonchem"]["version"] != 10
        error("CommonChem version other than 10 is not supported")
    end
    if length(data["molecules"]) != 1
        error("Only single molecule data is supported")
    end
    mol = data["molecules"][1]
    if mol["extensions"][1]["name"] != "rdkitRepresentation"
        error("Invalid RDKit CommonChem file format")
    end
    if mol["extensions"][1]["formatVersion"] != 2
        error("Unsupported RDKit CommonChem version")
    end
    return molgraph_from_dict(
        Val(:rdkit), T, CommonChemAtom, CommonChemBond, mol; kwargs...)
end


function to_dict(fmt::Val{:rdkit}, mol::MolGraph)
    data = Dict{String,Any}(
        "commonchem" => Dict("version" => 10),
        "defaults" => Dict(
            "atom" => Dict(
                "z" => 6,
                "impHs" => 0,
                "chg" => 0,
                "nRad" => 0,
                "isotope" => 0,
                "stereo" => "unspecified"
            ),
            "bond" => Dict(
                "bo" => 1,  # `type`` in CommonChem specification
                "stereo" => "unspecified"
            )
        ),
        "molecules" => [Dict(
            "atoms" => [],
            "bonds" => [],
            "extensions" => [Dict(
                "name" => "rdkitRepresentation",
                "formatVersion" => 2,
                "toolkitVersion" => "2022.09.5"
            )]
        )]
    )
    atomnum = atom_number(mol)
    implh = implicit_hydrogens(mol)
    chg = atom_charge(mol)
    mul = multiplicity(mol)
    ms = [atom_mass(props(mol, i)) for i in vertices(mol)]
    for i in vertices(mol)
        rcd = Dict{String,Any}()
        atomnum[i] == 6 || (rcd["z"] = atomnum[i])
        implh[i] == 0 || (rcd["impHs"] = implh[i])
        chg[i] == 0 || (rcd["chg"] = chg[i])
        mul[i] == 1 || (rcd["nRad"] = mul[i] - 1)
        isnothing(ms[i]) || (rcd["isotope"] = ms[i])
        if haskey(mol.gprops.stereocenter, i)
            center = mol.gprops.stereocenter[i]
            nbrs = ordered_neighbors(mol, i)
            rcd["stereo"] = isclockwise(center, nbrs[1:3]...) ? "cw" : "ccw"
        end
        push!(data["molecules"][1]["atoms"], rcd)
    end
    bondorder = bond_order(mol)
    for (i, e) in enumerate(edges(mol))
        rcd = Dict{String,Any}(
            "atoms" => [src(e) - 1, dst(e) - 1]
        )
        bondorder[i] == 1 || (rcd["bo"] = bondorder[i])
        if haskey(mol.gprops.stereobond, e)
            v1, v2, is_cis = mol.gprops.stereobond[e]
            rcd["stereoAtoms"] = [v1 - 1, v2 - 1]
            rcd["stereo"] = is_cis ? "cis" : "trans"
        end
        push!(data["molecules"][1]["bonds"], rcd)
    end
    data["molecules"][1]["conformers"] = []
    # Wedges cannot be preserved. Use coordgen to generate 2D coords manually
    for cds in mol.gprops.coords3d
        push!(
            data["molecules"][1]["conformers"],
            Dict("dim" => 3, "coords" => to_dict(Val(:default), Val(:coords3d), mol.gprops)))
    end
    return data
end


"""
    to_rdkdict(mol::MolGraph) -> Dict{String,Any}

Convert the molecule object into `Dict` in RDKit compatible CommonChem JSON format.
"""
to_rdkdict(x) = to_dict(Val{:rdkit}(), x)


"""
    to_rdkjson(mol::MolGraph) -> String

Convert the molecule object into RDKit compatible CommonChem JSON text.
"""
to_rdkjson(x) = to_json(Val{:rdkit}(), x)


"""
    rdkitmol(mol::MolGraph) -> RDKitMinimalLib.Mol

Convert the molecule object into a RDKit molecule object that can be used in RDKitMinimalLib.jl
"""
function rdkitmol(mol::MolGraph)
    return get_mol(to_json(Val{:rdkit}(), mol))
end


"""
    smiles(mol::MolGraph) -> String

Return a SMILES string of the molecule object.
"""
function smiles(mol::MolGraph, details=nothing)
    return get_smiles(rdkitmol(mol), details)
end


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
    morgan_fp_vector(mol::RDKitMinimalLib.Mol; kwargs...) -> BitArray
    morgan_fp_vector(mol::MolGraph; kwargs...) -> BitArray

Return a Morgan fingerprint bit array
"""
morgan_fp_vector(mol::Mol; kwargs...
    ) = uint8vec_to_bitarray(get_morgan_fp_as_bytes(mol::Mol, kwargs...))
morgan_fp_vector(mol::MolGraph; kwargs...
    ) = morgan_fp_vector(rdkitmol(mol); kwargs...)

"""
    rdkit_fp_vector(mol::RDKitMinimalLib.Mol; kwargs...) -> BitArray
    rdkit_fp_vector(mol::MolGraph; kwargs...) -> BitArray

Return a RDKit fingerprint bit array
"""
rdkit_fp_vector(mol::Mol; kwargs...
    ) = uint8vec_to_bitarray(get_rdkit_fp_as_bytes(mol::Mol, kwargs...))
rdkit_fp_vector(mol::MolGraph; kwargs...
    ) = rdkit_fp_vector(rdkitmol(mol); kwargs...)

"""
    pattern_fp_vector(mol::RDKitMinimalLib.Mol; kwargs...) -> BitArray
    pattern_fp_vector(mol::MolGraph; kwargs...) -> BitArray

Return a pattern fingerprint bit array, a topological fingerprint
optimized for substructure screening
"""
pattern_fp_vector(mol::Mol; kwargs...
    ) = uint8vec_to_bitarray(get_pattern_fp_as_bytes(mol::Mol, kwargs...))
pattern_fp_vector(mol::MolGraph; kwargs...
    ) = pattern_fp_vector(rdkitmol(mol); kwargs...)

"""
    atom_pair_fp_vector(mol::RDKitMinimalLib.Mol; kwargs...) -> BitArray
    atom_pair_fp_vector(mol::MolGraph; kwargs...) -> BitArray

Return a atom pairs fingerprint bit array
"""
atom_pair_fp_vector(mol::Mol; kwargs...
    ) = uint8vec_to_bitarray(get_atom_pair_fp_as_bytes(mol::Mol, kwargs...))
atom_pair_fp_vector(mol::MolGraph; kwargs...
    ) = atom_pair_fp_vector(rdkitmol(mol); kwargs...)

"""
    topological_torsion_fp_vector(mol::RDKitMinimalLib.Mol; kwargs...) -> BitArray
    topological_torsion_fp_vector(mol::MolGraph; kwargs...) -> BitArray

Return a topological torsions fingerprint bit array
"""
topological_torsion_fp_vector(mol::Mol; kwargs...
    ) = uint8vec_to_bitarray(get_topological_torsion_fp_as_bytes(mol::Mol, kwargs...))
topological_torsion_fp_vector(mol::MolGraph; kwargs...
    ) = topological_torsion_fp_vector(rdkitmol(mol); kwargs...)


jaccard_index(a::BitVector, b::BitVector
    ) = sum(a .& b) / sum(a .| b)  # Not exported