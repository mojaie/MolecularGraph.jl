#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

function rdk_on_init!(mol::SimpleMolGraph)
    # Do nothing
end


function rdk_on_update!(mol::SimpleMolGraph)
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


function reactive_molgraph(
        ::Type{T}, ::Type{V}, ::Type{E}, ::Type{G},
        data::JSON.Object{String,Any}) where {T,V<:CommonChemAtom,E<:CommonChemBond,G}
    gps = G()
    # edges
    es = Edge{T}[]
    eps = Dict{EdgeKey{T},E}()
    iscis = Dict("cis" => true, "trans" => false)
    for b in data["bonds"]
        edge = u_edge(T, b["atoms"][1] + 1, b["atoms"][2] + 1)
        push!(es, edge)
        if haskey(b, "stereo") && b["stereo"] in keys(iscis)
            u, v = b["stereoAtoms"]
            gps.stereobond[edge] = Stereobond{T}(u + 1, v + 1, iscis[b["stereo"]])
        end
        delete!(b, "atoms")
        delete!(b, "stereo")
        delete!(b, "stereoAtoms")
        eps[edge] = E(b)
    end
    g = SimpleGraph(es)
    # vertices
    vps = Dict{VertexKey{T},V}()
    isclockwise = Dict("cw" => true, "ccw" => false)
    for (i, a) in enumerate(data["atoms"])
        vps[i] = V(a)
        if haskey(a, "stereo") && a["stereo"] in keys(isclockwise)
            nbrs = ordered_neighbors(g, i)
            # TODO: ambiguity in implicit H
            gps.stereocenter[i] = Stereocenter{T}(nbrs[1:3]..., isclockwise[a["stereo"]])
        end
    end
    # expand fadjlist for vprops of isolated nodes
    for _ in nv(g):(length(vps) - 1)
        push!(g.fadjlist, T[])
    end
    if haskey(data, "conformers")
        for cds in data["conformers"]
            # Wedges cannot be preserved. Use coordgen to generate 2D coords manually
            cds["dim"] == 3 || continue
            push!(gps.descriptors.coords3d, [Point3d(cd...) for cd in cds["coords"]])
        end
    end
    return reactive_molgraph(g, vps, eps, gps)
end


function StructUtils.lift(::Type{MolGraph{T,CommonChemAtom,CommonChemBond}}, x) where T
    if x["commonchem"]["version"] != 10
        error("CommonChem version other than 10 is not supported")
    end
    if length(x["molecules"]) != 1
        error("Only single molecule data is supported")
    end
    moldata = x["molecules"][1]
    if moldata["extensions"][1]["name"] != "rdkitRepresentation"
        error("Invalid RDKit CommonChem file format")
    end
    if moldata["extensions"][1]["formatVersion"] != 2
        error("Unsupported RDKit CommonChem version")
    end
    mol = MolGraph{T,CommonChemAtom,CommonChemBond}(
        reactive_molgraph(T, CommonChemAtom, CommonChemBond, MolProperty{T}, x)...)
    mol.state.on_init = rdk_on_init!
    mol.state.on_update = rdk_on_update!
    return mol
end


"""
    to_rdkdict(mol::MolGraph) -> Dict{String,Any}

Convert the molecule object into `Dict` in RDKit compatible CommonChem JSON format.
"""
function to_rdkdict(mol::SimpleMolGraph)
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
    iso = isotope(mol)
    for i in vertices(mol)
        rcd = Dict{String,Any}()
        atomnum[i] == 6 || (rcd["z"] = atomnum[i])
        implh[i] == 0 || (rcd["impHs"] = implh[i])
        chg[i] == 0 || (rcd["chg"] = chg[i])
        mul[i] == 1 || (rcd["nRad"] = mul[i] - 1)
        iso[i] == 0 || (rcd["isotope"] = iso[i])
        if haskey(mol[:stereocenter], i)
            center = mol[:stereocenter][i]
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
        if haskey(mol[:stereobond], e)
            v1, v2, is_cis = mol[:stereobond][e]
            rcd["stereoAtoms"] = [v1 - 1, v2 - 1]
            rcd["stereo"] = is_cis ? "cis" : "trans"
        end
        push!(data["molecules"][1]["bonds"], rcd)
    end
    data["molecules"][1]["conformers"] = []
    # Wedges cannot be preserved. Use coordgen to generate 2D coords manually
    for crds in mol[:descriptors].coords3d
        push!(
            data["molecules"][1]["conformers"],
            Dict("dim" => 3, "coords" => [[p...] for p in crds]))
    end
    return data
end


"""
    to_rdkjson(mol::MolGraph) -> String

Convert the molecule object into RDKit compatible CommonChem JSON text.
"""
to_rdkjson(x) = JSON.json(x)


"""
    rdkitmol(mol::MolGraph) -> RDKitMinimalLib.Mol

Convert the molecule object into a RDKit molecule object that can be used in RDKitMinimalLib.jl
"""
function rdkitmol end


"""
    smiles(mol::MolGraph) -> String

Return a SMILES string of the molecule object.
"""
function smiles end


"""
    morgan_fp_vector(mol::MolGraph, details=nothing) -> BitArray
    morgan_fp_string(mol::MolGraph, details=nothing) -> String

Return a Morgan fingerprint bit array
"""
function morgan_fp_vector end
function morgan_fp_string end


"""
    rdkit_fp_vector(mol::MolGraph, details=nothing) -> BitArray
    rdkit_fp_string(mol::MolGraph, details=nothing) -> String

Return a RDKit fingerprint bit array
"""
function rdkit_fp_vector end
function rdkit_fp_string end


"""
    pattern_fp_vector(mol::MolGraph, details=nothing) -> BitArray
    pattern_fp_string(mol::MolGraph, details=nothing) -> String

Return a pattern fingerprint bit array, a topological fingerprint
optimized for substructure screening
"""
function pattern_fp_vector end
function pattern_fp_string end


"""
    atom_pair_fp_vector(mol::MolGraph, details=nothing) -> BitArray
    atom_pair_fp_string(mol::MolGraph, details=nothing) -> String

Return a atom pairs fingerprint bit array
"""
function atom_pair_fp_vector end
function atom_pair_fp_string end


"""
    topological_torsion_fp_vector(mol::MolGraph, details=nothing) -> BitArray
    topological_torsion_fp_string(mol::MolGraph, details=nothing) -> String

Return a topological torsions fingerprint bit array
"""
function topological_torsion_fp_vector end
function topological_torsion_fp_string end



jaccard_index(a::BitVector, b::BitVector
    ) = sum(a .& b) / sum(a .| b)  # Not exported
