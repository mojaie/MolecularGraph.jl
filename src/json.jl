#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    to_dict, to_json


"""
    to_dict(mol::MolGraph) -> Dict{String,Any}

Convert molecule object into JSON compatible dictionary.
"""
to_dict(::Val, mol::AbstractMolGraph) = error("method to_dict not implemented")
to_dict(::Val, atom_or_bond::Dict) = atom_or_bond
to_dict(::Val, metadata::AbstractString) = metadata
to_dict(::Val, value::Number) = value
to_dict(::Val{:default}, meta::Metadata) = [[i, val] for (i, val) in meta]
to_dict(x) = to_dict(Val{:default}(), x)


function to_dict(fmt::Val{:default}, mol::MolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    return Dict(
        "eltype" => string(eltype(mol)),
        "vproptype" => string(vproptype(mol)),
        "eproptype" => string(eproptype(mol)),
        "graph" => [[src(e), dst(e)] for e in edges(mol)],
        "vprops" => [to_dict(fmt, props(mol, i)) for i in vertices(mol)],
        "eprops" => [to_dict(fmt, props(mol, e)) for e in edges(mol)],
        "gprops" => [[string(k), string(typeof(v)), applicable(to_dict, fmt, v) ? to_dict(fmt, v) : v] for (k, v) in mol.gprops],
        "caches" => [[string(k), string(typeof(v)), applicable(to_dict, fmt, v) ? to_dict(fmt, v) : v] for (k, v) in mol.state[:caches]]
    )
end


"""
    to_json(mol::MolGraph) -> String

Convert molecule object into JSON String.
"""
to_json(fmt::Val, mol::AbstractMolGraph) = JSON.json(to_dict(fmt, mol))
to_json(x) = to_json(Val{:default}(), x)


function molgraph_from_dict(
            ::Val{:default}, ::Type{T}, data::Dict, config::Dict{Symbol,Any}; kwargs...
        ) where T <: AbstractMolGraph
    I = eltype(T)
    V = vproptype(T)
    E = eproptype(T)
    g = SimpleGraph(Edge{I}[Edge{I}(e...) for e in data["graph"]])
    vps = Dict{I,V}(i => V(vp) for (i, vp) in enumerate(data["vprops"]))
    eps = Dict{Edge{I},E}(e => E(ep) for (e, ep) in zip(edges(g), data["eprops"]))
    gps = Dict{Symbol,Any}(
        Symbol(key) => eval(Meta.parse(dtype))(prop) for (key, dtype, prop) in data["gprops"])
    default_config = Dict{Symbol,Any}(
        :caches => Dict{Symbol,Any}(
            Symbol(key) => eval(Meta.parse(dtype))(prop) for (key, dtype, prop) in data["caches"]),
        :initialized => true,
        :has_updates => false
    )
    merge!(default_config, config)
    mol = T(g, vps, eps, gprop_map=gps, config_map=default_config; kwargs...)
    update_edge_rank!(mol)  # initialization skiped, but edge_rank should be resumed
    return mol
end


# JSON auto detect

function MolGraph(data::Dict; config=Dict{Symbol,Any}(), kwargs...)
    default_config = Dict{Symbol,Any}()
    if haskey(data, "commonchem")
        T = Int
        V = CommonChemAtom
        E = CommonChemBond
        fmt = Val{:rdkit}()
        default_config = Dict{Symbol,Any}(:on_init => rdk_on_init!, :updater => rdk_on_update!)
    elseif haskey(data, "eltype")
        T = eval(Meta.parse(data["eltype"]))
        V = eval(Meta.parse(data["vproptype"]))
        E = eval(Meta.parse(data["eproptype"]))
        fmt = Val{:default}()
        if V == SDFAtom && E == SDFBond
            default_config = Dict{Symbol,Any}(:on_init => sdf_on_init!, :updater => sdf_on_update!)
        elseif V == SMILESAtom && E == SMILESBond
            default_config = Dict{Symbol,Any}(:on_init => smiles_on_init!, :updater => smiles_on_update!)
        end
    else
        error("Invalid JSON format")
    end
    merge!(default_config, config)
    return molgraph_from_dict(fmt, MolGraph{T,V,E}, data, default_config; kwargs...)
end

MolGraph(json::String; kwargs...) = MolGraph(JSON.parse(json); kwargs...)




