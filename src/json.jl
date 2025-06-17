#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

"""
    to_dict(fmt::Val{:default}, mol::MolGraph) -> Dict{String,Any}

Convert molecule object into JSON compatible dictionary.
"""
function to_dict(fmt::Val{:default}, mol::ReactiveMolGraph{T,V,E}) where {T,V,E}
    dispatch_update!(mol)
    return Dict(
        "vproptype" => string(nameof(V)),
        "eproptype" => string(nameof(E)),
        "graph" => [[src(e), dst(e)] for e in edges(mol)],
        "vprops" => [to_dict(fmt, props(mol, i)) for i in vertices(mol)],
        "eprops" => [to_dict(fmt, props(mol, e)) for e in edges(mol)],
        "gprops" => to_dict(fmt, mol.gprops)
    )
end

to_dict(mol::AbstractMolGraph) = to_dict(Val{:default}(), mol)

"""
    to_json(fmt::Val, mol::AbstractMolGraph) -> String
    to_json(mol::AbstractMolGraph) -> String

Convert molecule object into JSON String.
"""
to_json(fmt::Val, mol::AbstractMolGraph) = JSON.json(to_dict(fmt, mol))
to_json(mol::AbstractMolGraph) = to_json(Val{:default}(), mol)


function reactive_molgraph(
        ::Val{:default}, ::Type{T}, ::Type{V}, ::Type{E},
        data::Dict, config::MolState) where {T,V,E}
    g = SimpleGraph(Edge{T}[Edge{T}(e...) for e in data["graph"]])
    vps = Dict{T,V}(i => V(vp) for (i, vp) in enumerate(data["vprops"]))
    eps = Dict{Edge{T},E}(e => E(ep) for (e, ep) in zip(edges(g), data["eprops"]))
    gps = MolProperty{T}(data["gprops"])
    return (g, vps, eps, gps, config)
end


function MolGraph{T,SDFAtom,SDFBond}(data::Dict
        ; on_init=sdf_on_init!, on_update=sdf_on_update!) where T<:Integer
    if data["vproptype"] != "SDFAtom" || data["eproptype"] != "SDFBond"
        error("Incompatible element property types")
    end
    config=MolState{T}(;
        on_init=on_init,
        on_update=on_update,
        initialized = true,  # Skip initialization
        has_updates = false  # Do not update ready-to-use descriptors!
    )
    mol = MolGraph(
        reactive_molgraph(Val(:default), T, SDFAtom, SDFBond, data, config)...)
    update_edge_rank!(mol)  # but edge_rank should be resumed
    return mol
end

function MolGraph{T,SMILESAtom,SMILESBond}(data::Dict
        ; on_init=smiles_on_init!, on_update=smiles_on_update!) where T<:Integer
    if data["vproptype"] != "SMILESAtom" || data["eproptype"] != "SMILESBond"
        error("Incompatible element property types")
    end
    config=MolState{T}(;
        on_init=on_init,
        on_update=on_update,
        initialized = true,  # Skip initialization
        has_updates = false  # Do not update ready-to-use descriptors!
    )
    mol = MolGraph(
        reactive_molgraph(Val(:default), T, SMILESAtom, SMILESBond, data, config)...)
    update_edge_rank!(mol)  # but edge_rank should be resumed
    return mol
end

function QueryMolGraph{T,QueryAtom,QueryBond}(data::Dict
        ; on_init=default_on_init!, on_update=default_on_update!) where T<:Integer
    if data["vproptype"] != "QueryAtom" || data["eproptype"] != "QueryBond"
        error("Incompatible element property types")
    end
    config=MolState{T}(
        on_init=on_init,
        on_update=on_update,
        initialized = true,  # Skip initialization
        has_updates = false  # Do not update ready-to-use descriptors!
    )
    mol = QueryMolGraph(
        reactive_molgraph(Val(:default), T, QueryAtom, QueryBond, data, config)...)
    update_edge_rank!(mol)  # but edge_rank should be resumed
    return mol
end


# JSON auto detect

function MolGraph(data::Dict; kwargs...)
    if haskey(data, "commonchem")
        return MolGraph{Int,CommonChemAtom,CommonChemBond}(data; kwargs...)
    end
    if !haskey(data, "vproptype") || !haskey(data, "eproptype")
        error("Invalid JSON format")
    end
    if data["vproptype"] == "SDFAtom" && data["eproptype"] == "SDFBond"
        return MolGraph{Int,SDFAtom,SDFBond}(data; kwargs...)
    elseif data["vproptype"] == "SMILESAtom" && data["eproptype"] == "SMILESBond"
        return MolGraph{Int,SMILESAtom,SMILESBond}(data; kwargs...)
    end
    error("Invalid JSON format")
end

function QueryMolGraph(data::Dict; kwargs...)
    if !haskey(data, "vproptype") || !haskey(data, "eproptype")
        error("Invalid JSON format")
    end
    if data["vproptype"] == "QueryAtom" && data["eproptype"] == "QueryBond"
        return QueryMolGraph{Int,QueryAtom,QueryBond}(data; kwargs...)
    end
    error("Invalid JSON format")
end


MolGraph(json::String; kwargs...) = MolGraph(JSON.parse(json); kwargs...)
QueryMolGraph(json::String; kwargs...) = QueryMolGraph(JSON.parse(json); kwargs...)

MolGraph{T,V,E}(json::String; kwargs...
    ) where {T,V,E} = MolGraph{T,V,E}(JSON.parse(json); kwargs...)
QueryMolGraph{T,V,E}(json::String; kwargs...
    ) where {T,V,E} = QueryMolGraph{T,V,E}(JSON.parse(json); kwargs...)
