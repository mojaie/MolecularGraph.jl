#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    to_dict, to_json, to_rdkdict, to_rdkjson, to_rdkmol, to_rdkqmol


"""
    to_dict(mol::MolGraph) -> Dict{String,Any}

Convert molecule object into JSON compatible dictionary.
"""
to_dict(::Val, mol::AbstractMolGraph) = error("method to_dict not implemented")
to_dict(::Val, atom_or_bond::Dict) = atom_or_bond
to_dict(::Val, metadata::AbstractString) = metadata
to_dict(::Val, value::Number) = value
to_dict(x) = to_dict(Val{:standard}(), x)
to_rdkdict(x) = to_dict(Val{:rdkit}(), x)


"""
    to_json(mol::MolGraph) -> String

Convert molecule object into JSON String.
"""
to_json(fmt::Val, mol::AbstractMolGraph) = JSON.json(to_dict(fmt, mol))
to_json(x) = to_json(Val{:standard}(), x)
to_rdkjson(x) = to_json(Val{:rdkit}(), x)


function to_dict(::Val{:rdkit}, mol::MolGraph)
    jdict = Dict{String,Any}(
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
                "bo" => 1,
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
        rcd = Dict{String,Any}()
        get_prop(mol, i, :symbol) === :C || (rcd["z"] = atomnumber(get_prop(mol, i, :symbol)))
        implh_[i] == 0 || (rcd["impHs"] = implh_[i])
        get_prop(mol, i, :charge) == 0 || (rcd["chg"] = get_prop(mol, i, :charge))
        get_prop(mol, i, :multiplicity) == 1 || (rcd["nRad"] = get_prop(mol, i, :multiplicity) - 1)
        isnothing(get_prop(mol, i, :mass)) || (rcd["isotope"] = get_prop(mol, i, :mass))
        get_prop(mol, i, :stereo) === :unspecified || (rcd["stereo"] = get_prop(mol, i, :stereo))
        push!(jdict["molecules"][1]["atoms"], rcd)
    end
    for (i, e) in enumerate(edges(mol))
        rcd = Dict{String,Any}("atoms" => [])
        get_prop(mol, e, :order) == 1 || (rcd["bo"] = get_prop(mol, e, :order))
        get_prop(mol, i, :stereo) === :unspecified || (rcd["stereo"] = get_prop(mol, i, :stereo))
        push!(jdict["molecules"][1]["bonds"], rcd)
    end
    return jdict
end


function molgraph_from_rdkdict(
        T::Type{<:Integer}, V::Type, E::Type, data::Dict, config::Dict{Symbol,Any})

end

function to_rdkmol(mol::MolGraph)
    return get_mol(to_json(Val{:rdkit}(), mol))
end


function to_rdkqmol(mol::MolGraph{Int,QueryTree,QueryTree})
    return get_qmol(to_json(Val{:rdkit}(), mol))
end
