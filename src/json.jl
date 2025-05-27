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
        "molecules" => [
            Dict(
                "atoms" => [],
                "bonds" => [],
                "extensions" => [
                    Dict(
                        "name" => "rdkitRepresentation",
                        "formatVersion" => 2,
                        "toolkitVersion" => "2022.09.5"
                    )
                ]
            )
        ]
    )
    implh_ = implicit_hydrogens(mol)
    for i in vertices(mol)
        rcd = to_dict(fmt, props(mol, i))
        implh_[i] == 0 || (rcd["impHs"] = implh_[i])
        # get_prop(mol, i, :stereo) === :unspecified || (rcd["stereo"] = get_prop(mol, i, :stereo))
        push!(data["molecules"][1]["atoms"], rcd)
    end
    for (i, e) in enumerate(edges(mol))
        rcd = to_dict(fmt, props(mol, e))
        rcd["atoms"] = [src(e) - 1, dst(e) - 1]
        # get_prop(mol, i, :stereo) === :unspecified || (rcd["stereo"] = get_prop(mol, i, :stereo))
        push!(data["molecules"][1]["bonds"], rcd)
    end
    return data
end


function rdktomol(::Type{T}, data::Dict) where T <: AbstractMolGraph
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
    # coords
    # stereo
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