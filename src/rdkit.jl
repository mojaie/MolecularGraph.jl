#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    to_rdkdict, to_rdkjson,
    rdkitmol, smiles,
    morgan_fp_vector, rdkit_fp_vector,
    pattern_fp_vector, atom_pair_fp_vector,
    topological_torsion_fp_vector

using RDKitMinimalLib: Mol, get_mol, get_smiles,
    get_morgan_fp_as_bytes, get_rdkit_fp_as_bytes,
    get_pattern_fp_as_bytes, get_atom_pair_fp_as_bytes,
    get_topological_torsion_fp_as_bytes


function rdk_on_init!(mol)
    update_edge_rank!(mol)
    set_state!(mol, :initialized, true)
end


function rdk_on_update!(mol)
    update_edge_rank!(mol)
    clear_caches!(mol)
    set_state!(mol, :has_updates, false)
    # preprocessing
    #  (add some preprocessing methods here)
    # recalculate bottleneck descriptors
    sssr!(mol)
    lone_pair!(mol)
    apparent_valence!(mol)
    valence!(mol)
    is_ring_aromatic!(mol)
end


function molgraph_from_dict(
            ::Val{:rdkit}, ::Type{T}, data::Dict, config::Dict{Symbol,Any}; kwargs...
        ) where T <: AbstractMolGraph
    data["commonchem"]["version"] == 10 || error("CommonChem version other than 10 is not supported")
    mol = data["molecules"][1]  # only single mol data is supported
    mol["extensions"][1]["name"] == "rdkitRepresentation" || error("Invalid RDKit CommonChem file format")
    mol["extensions"][1]["formatVersion"] == 2 || error("Unsupported RDKit CommonChem version")
    I = eltype(T)
    V = vproptype(T)
    E = eproptype(T)
    gps = Dict{Symbol,Any}(
        :metadata => Metadata()
    )
    # edges
    es = Edge{I}[]
    eps = Dict{Edge{I},E}()
    stereobonds = Dict{Edge{I},Tuple{I,I,Bool}}()
    iscis = Dict("cis" => true, "trans" => false)
    for b in mol["bonds"]
        edge = u_edge(I, b["atoms"][1] + 1, b["atoms"][2] + 1)
        push!(es, edge)
        if haskey(b, "stereo") && b["stereo"] in keys(iscis)
            u, v = b["stereoAtoms"]
            stereobonds[edge] = (u + 1, v + 1, iscis[b["stereo"]])
        end
        eps[edge] = CommonChemBond(b)
    end
    g = SimpleGraph(es)
    gps[:stereobond] = Stereobond{I}(stereobonds)
    # vertices
    vps = Dict{I,V}()
    stereocenters = Dict{I,Tuple{I,I,I,Bool}}()
    isclockwise = Dict("cw" => true, "ccw" => false)
    for (i, a) in enumerate(mol["atoms"])
        vps[i] = CommonChemAtom(a)
        if haskey(a, "stereo") && a["stereo"] in keys(isclockwise)
            nbrs = ordered_neighbors(g, i)
            # TODO: ambiguity in implicit H
            stereocenters[i] = (nbrs[1:3]..., isclockwise[a["stereo"]])
        end
    end
    gps[:stereocenter] = Stereocenter{I}(stereocenters)
    # TODO: multiple coords
    # TODO: wedge bond notation for drawing not supported yet, use coordgen
    if haskey(mol, "conformers")
        coords = zeros(Float64, nv(g), 2)
        for (i, c) in enumerate(mol["conformers"][1]["coords"])
            coords[i, :] = c[1:2]
        end
        gps[:coords2d] = coords
    end
    return T(g, vps, eps, gprop_map=gps, config_map=config)
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
    chg = charge(mol)
    mul = multiplicity(mol)
    ms = [mass(props(mol, i)) for i in vertices(mol)]
    stereocenters = has_prop(mol, :stereocenter) ? get_prop(mol, :stereocenter) : Dict()
    for i in vertices(mol)
        rcd = Dict{String,Any}()
        atomnum[i] == 6 || (rcd["z"] = atomnum[i])
        implh[i] == 0 || (rcd["impHs"] = implh[i])
        chg[i] == 0 || (rcd["chg"] = chg[i])
        mul[i] == 1 || (rcd["nRad"] = mul[i] - 1)
        isnothing(ms[i]) || (rcd["isotope"] = ms[i])
        if haskey(stereocenters, i)
            center = stereocenters[i]
            nbrs = ordered_neighbors(mol, i)
            rcd["stereo"] = isclockwise(center, nbrs[1:3]...) ? "cw" : "ccw"
        end
        push!(data["molecules"][1]["atoms"], rcd)
    end
    stereobonds = has_prop(mol, :stereobond) ? get_prop(mol, :stereobond) : Dict()
    bondorder = bond_order(mol)
    for (i, e) in enumerate(edges(mol))
        rcd = Dict{String,Any}(
            "atoms" => [src(e) - 1, dst(e) - 1]
        )
        bondorder[i] == 1 || (rcd["bo"] = bondorder[i])
        if haskey(stereobonds, e)
            v1, v2, is_cis = stereobonds[e]
            rcd["stereoAtoms"] = [v1 - 1, v2 - 1]
            rcd["stereo"] = is_cis ? "cis" : "trans"
        end
        push!(data["molecules"][1]["bonds"], rcd)
    end
    # TODO: multiple coords
    if has_coords(mol)
        data["molecules"][1]["conformers"] = []
        coords_ = coords2d(mol)
        push!(
            data["molecules"][1]["conformers"],
            Dict("dim" => 2, "coords" => [coords_[i, 1:2] for i in vertices(mol)])
        )
    end
    return data
end


to_rdkdict(x) = to_dict(Val{:rdkit}(), x)
to_rdkjson(x) = to_json(Val{:rdkit}(), x)


function rdkitmol(mol::MolGraph)
    return get_mol(to_json(Val{:rdkit}(), mol))
end


function smiles(mol::MolGraph, details=nothing)
    return get_smiles(rdkitmol(mol), details)
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
    ) = morgan_fp_vector(rdkitmol(mol); kwargs...)

rdkit_fp_vector(mol::Mol; kwargs...
    ) = uint8vec_to_bitarray(get_rdkit_fp_as_bytes(mol::Mol, kwargs...))
rdkit_fp_vector(mol::MolGraph; kwargs...
    ) = rdkit_fp_vector(rdkitmol(mol); kwargs...)

pattern_fp_vector(mol::Mol; kwargs...
    ) = uint8vec_to_bitarray(get_pattern_fp_as_bytes(mol::Mol, kwargs...))
pattern_fp_vector(mol::MolGraph; kwargs...
    ) = pattern_fp_vector(rdkitmol(mol); kwargs...)

atom_pair_fp_vector(mol::Mol; kwargs...
    ) = uint8vec_to_bitarray(get_atom_pair_fp_as_bytes(mol::Mol, kwargs...))
atom_pair_fp_vector(mol::MolGraph; kwargs...
    ) = atom_pair_fp_vector(rdkitmol(mol); kwargs...)

topological_torsion_fp_vector(mol::Mol; kwargs...
    ) = uint8vec_to_bitarray(get_topological_torsion_fp_as_bytes(mol::Mol, kwargs...))
topological_torsion_fp_vector(mol::MolGraph; kwargs...
    ) = topological_torsion_fp_vector(rdkitmol(mol); kwargs...)


jaccard_index(a::BitVector, b::BitVector
    ) = sum(a .& b) / sum(a.| b)  # Not exported