#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    to_dict, to_json, to_rdkdict, to_rdkjson, rdktomol


const CommonChemMolGraph = MolGraph{Int,CommonChemAtom,CommonChemBond}


"""
    to_dict(mol::MolGraph) -> Dict{String,Any}

Convert molecule object into JSON compatible dictionary.
"""
to_dict(::Val, mol::AbstractMolGraph) = error("method to_dict not implemented")
to_dict(::Val, atom_or_bond::Dict) = atom_or_bond
to_dict(::Val, metadata::AbstractString) = metadata
to_dict(::Val, value::Number) = value
to_dict(x) = to_dict(Val{:default}(), x)
to_rdkdict(x) = to_dict(Val{:rdkit}(), x)


"""
    to_json(mol::MolGraph) -> String

Convert molecule object into JSON String.
"""
to_json(fmt::Val, mol::AbstractMolGraph) = JSON.json(to_dict(fmt, mol))
to_json(x) = to_json(Val{:default}(), x)
to_rdkjson(x) = to_json(Val{:rdkit}(), x)


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
    stereocenters = get_prop(mol, :stereocenter)
    for i in vertices(mol)
        rcd = Dict{String,Any}()
        atomnum[i] == 6 || (rcd["z"] = atomnum[i])
        implh[i] == 0 || (rcd["impHs"] = implh[i])
        chg[i] == 0 || (rcd["chg"] = chg[i])
        mul[i] == 0 || (rcd["nRad"] = mul[i] - 1)
        isnothing(ms[i]) || (rcd["isotope"] = ms[i])
        if haskey(stereocenters, i)
            center = stereocenters[i]
            nbrs = ordered_neighbors(mol, i)
            rcd["stereo"] = isclockwise(center, nbrs[1:3]...) ? "cw" : "ccw"
        end
        push!(data["molecules"][1]["atoms"], rcd)
    end
    stereobonds = get_prop(mol, :stereobond)
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
    # TODO: coords 2d and 3d
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


function rdktomol(::Type{T}, data::Dict) where T <: AbstractMolGraph
    error("not implemented yet")
    data["commonchem"]["version"] == 10 || error("CommonChem version other than 10 is not supported")
    mol = data["molecules"][1]  # only single mol data is supported
    mol["extensions"][1]["name"] == "rdkitRepresentation" || error("Invalid RDKit CommonChem file format")
    mol["extensions"][1]["formatVersion"] == 2 || error("Unsupported RDKit CommonChem version")
    I = eltype(T)
    V = vproptype(T)
    E = eproptype(T)
    es = Edge{I}[]
    vps = Dict{I,V}()
    eps = Dict{Edge{I},E}()
    # TODO: coords
    # TODO: stereo
    gps = Dict()
    for (i, a) in enumerate(mol["atoms"])
        vps[i] = CommonChemAtom(a)
    end
    for b in mol["bonds"]
        edge = u_edge(I, b["atoms"][1] + 1, b["atoms"][2] + 1)
        print(edge)
        push!(es, edge)
        eps[edge] = CommonChemBond(b)
    end
    return MolGraph{I,V,E}(SimpleGraph(es), vps, eps, gprop_map=gps)
end


rdktomol(data::Dict) = rdktomol(CommonChemMolGraph, data)
rdktomol(json::String) = rdktomol(CommonChemMolGraph, JSON.parse(json))
