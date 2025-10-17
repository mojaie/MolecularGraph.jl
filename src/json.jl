#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

function JSON.lower(x::ReactiveMolGraph{T,V,E}) where {T,V,E}
    dispatch_update!(x)
    return Dict{String,Any}(
        "vproptype" => string(nameof(V)),
        "eproptype" => string(nameof(E)),
        "graph" => x.graph,
        "vprops" => x.vprops,
        "eprops" => x.eprops,
        "gprops" => x.gprops,
        "state" => Dict{String,Any}()
    )
end


function MolGraph{T,V,E}(json::AbstractString
        ; on_init=sdf_on_init!, on_update=sdf_on_update!) where T<:Integer
    mol = JSON.parse(data, MolGraph{T,V,E})
    mol.config.on_init = on_init
    mol.config.on_update = on_update
    mol.config.initialized = true,  # Skip initialization
    mol.config.has_updates = false  # Do not update ready-to-use descriptors!
    return mol
end




function MolGraph{T,SMILESAtom,SMILESBond}(data::JSON.Object{String,Any}
        ; on_init=smiles_on_init!, on_update=smiles_on_update!) where T<:Integer
    if data["vproptype"] != "SMILESAtom" || data["eproptype"] != "SMILESBond"
        error("Incompatible element property types")
    end
    gps = reconstruct(MolProperty{T}, data["gprops"])
    config=MolState{T}(;
        on_init=on_init,
        on_update=on_update,
        initialized = true,  # Skip initialization
        has_updates = false  # Do not update ready-to-use descriptors!
    )
    mol = MolGraph(
        reactive_molgraph(Val(:default), T, SMILESAtom, SMILESBond, data, gps, config)...)
    return mol
end

function QueryMolGraph{T,QueryAtom,QueryBond}(data::JSON.Object{String,Any}
        ; on_init=default_on_init!, on_update=default_on_update!) where T<:Integer
    if data["vproptype"] != "QueryAtom" || data["eproptype"] != "QueryBond"
        error("Incompatible element property types")
    end
    gps = reconstruct(QueryMolProperty{T}, data["gprops"])
    config=MolState{T}(
        on_init=on_init,
        on_update=on_update,
        initialized = true,  # Skip initialization
        has_updates = false  # Do not update ready-to-use descriptors!
    )
    mol = QueryMolGraph(
        reactive_molgraph(Val(:default), T, QueryAtom, QueryBond, data, gps, config)...)
    return mol
end


# JSON auto detect

function mol_from_dict(data::JSON.Object{String,Any}; kwargs...)
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
    elseif data["vproptype"] == "QueryAtom" && data["eproptype"] == "QueryBond"
        return QueryMolGraph{Int,QueryAtom,QueryBond}(data; kwargs...)
    end
    error("Invalid JSON format")
end


function MolGraph(data::JSON.Object{String,Any}; kwargs...)
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

function QueryMolGraph(data::JSON.Object{String,Any}; kwargs...)
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
