#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

"""
    to_dict(fmt::Val{:default}, mol::MolGraph) -> Dict{String,Any}

Convert molecule object into JSON compatible dictionary.
"""
function to_dict(fmt::Val{:default}, mol::MolGraph)
    dispatch_update!(mol)
    rev_type = Dict(v => k for (k, v) in ELEMENT_TYPE_REGISTRY)
    return Dict(
        "eltype" => get(rev_type, eltype(mol), Int),
        "vproptype" => rev_type[vproptype(mol)],
        "eproptype" => rev_type[eproptype(mol)],
        "graph" => [[src(e), dst(e)] for e in edges(mol)],
        "vprops" => [to_dict(fmt, props(mol, i)) for i in vertices(mol)],
        "eprops" => [to_dict(fmt, props(mol, e)) for e in edges(mol)],
        "gprops" => to_dict(fmt, mol.gprops)
    )
end


"""
    to_json(fmt::Val, mol::AbstractMolGraph) -> String
    to_json(mol::AbstractMolGraph) -> String

Convert molecule object into JSON String.
"""
to_json(fmt::Val, mol::AbstractMolGraph) = JSON.json(to_dict(fmt, mol))
to_json(mol::AbstractMolGraph) = to_json(Val{:default}(), mol)


function molgraph_from_dict(
        ::Val{:default}, ::Type{T}, ::Type{V}, ::Type{E},
        data::Dict, config::MolGraphState) where {T,V,E}
    g = SimpleGraph(Edge{T}[Edge{T}(e...) for e in data["graph"]])
    vps = Dict{T,V}(i => V(vp) for (i, vp) in enumerate(data["vprops"]))
    eps = Dict{Edge{T},E}(e => E(ep) for (e, ep) in zip(edges(g), data["eprops"]))
    gps = MolGraphProperty{T}(data["gprops"])
    # Skip initialization to use stored descriptors without recalculation
    config.initialized = true
    config.has_updates = false
    mol = MolGraph{T,V,E}(g, vps, eps, gps, config)
    update_edge_rank!(mol)  # but edge_rank should be resumed
    return mol
end


function MolGraph{T,SDFAtom,SDFBond}(data::Dict;
        on_init=sdf_on_init!, on_update=sdf_on_update!) where T<:Integer
    if data["vproptype"] != "SDFAtom" || data["eproptype"] != "SDFBond"
        error("Incompatible element property types")
    end
    config=MolGraphState{T}(on_init, on_update)
    return molgraph_from_dict(
        Val(:default), T, SDFAtom, SDFBond, data, config)
end

function MolGraph{T,SMILESAtom,SMILESBond}(data::Dict;
        on_init=smiles_on_init!, on_update=smiles_on_update!) where T<:Integer
    if data["vproptype"] != "SMILESAtom" || data["eproptype"] != "SMILESBond"
        error("Incompatible element property types")
    end
    config=MolGraphState{T}(on_init, on_update)
    return molgraph_from_dict(
        Val(:default), T, SMILESAtom, SMILESBond, data, config)
end

function MolGraph{T,QueryTree,QueryTree}(data::Dict;
        on_init=default_on_init!, on_update=default_on_update!) where T<:Integer
    if data["vproptype"] != "QueryTree" || data["eproptype"] != "QueryTree"
        error("Incompatible element property types")
    end
    config=MolGraphState{T}(on_init, on_update)
    return molgraph_from_dict(
        Val(:default), T, QueryTree, QueryTree, data, config)
end


# JSON auto detect

function MolGraph(data::Dict; kwargs...)
    if haskey(data, "commonchem")
        return MolGraph{Int,CommonChemAtom,CommonChemBond}(data; kwargs...)
    end
    if !haskey(data, "eltype") || !haskey(data, "vproptype") || !haskey(data, "eproptype")
        error("Invalid JSON format")
    end
    if data["eltype"] == "Int"
        if data["vproptype"] == "SDFAtom" && data["eproptype"] == "SDFBond"
            return MolGraph{Int,SDFAtom,SDFBond}(data; kwargs...)
        elseif data["vproptype"] == "SMILESAtom" && data["eproptype"] == "SMILESBond"
            return MolGraph{Int,SMILESAtom,SMILESBond}(data; kwargs...)
        elseif data["vproptype"] == "QueryTree" && data["eproptype"] == "QueryTree"
            return MolGraph{Int,QueryTree,QueryTree}(data; kwargs...)
        end
    end
    error("Invalid JSON format")
end

MolGraph(json::String; kwargs...) = MolGraph(JSON.parse(json); kwargs...)
MolGraph{T,V,E}(json::String; kwargs...) where {T,V,E} = MolGraph{T,V,E}(JSON.parse(json); kwargs...)




