#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

function StructUtils.lower(x::ReactiveMolGraph{T,V,E}) where {T,V,E}
    dispatch_update!(x)  # update descriptors before serialization
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

StructUtils.structlike(::StructUtils.StructStyle, ::Type{MolGraph{T,V,E}}) where {T,V,E} = false
StructUtils.structlike(::StructUtils.StructStyle, ::Type{QueryMolGraph{T,V,E}}) where {T,V,E} = false

reactive_molgraph(
    ::Type{T}, ::Type{V}, ::Type{E}, ::Type{G}, x::JSON.Object{String,Any}
) where {T,V,E,G} = reactive_molgraph(
    StructUtils.make(SimpleGraph{T}, x["graph"]),
    StructUtils.make(Dict{VertexKey{T},V}, x["vprops"]),
    StructUtils.make(Dict{EdgeKey{T},E}, x["eprops"]),
    StructUtils.make(G, x["gprops"])
    ;initialized=true,  # Skip initialization
    has_updates=false  # Do not update ready-to-use descriptors!
)

function StructUtils.lift(::Type{MolGraph{T,SDFAtom,SDFBond}}, x) where T
    mol = MolGraph{T,SDFAtom,SDFBond}(reactive_molgraph(T, SDFAtom, SDFBond, MolProperty{T}, x)...)
    mol.state.on_init = sdf_on_init!
    mol.state.on_update = sdf_on_update!
    return mol
end

function StructUtils.lift(::Type{MolGraph{T,SMILESAtom,SMILESBond}}, x) where T
    mol = MolGraph{T,SMILESAtom,SMILESBond}(reactive_molgraph(T, SMILESAtom, SMILESBond, MolProperty{T}, x)...)
    mol.state.on_init = smiles_on_init!
    mol.state.on_update = smiles_on_update!
    return mol
end

StructUtils.lift(::Type{QueryMolGraph{T,V,E}}, x
    ) where {T,V,E} = QueryMolGraph{T,V,E}(reactive_molgraph(T, V, E, QueryMolProperty{T}, x)...)


function MolGraph{T,V,E}(json::AbstractString; on_init=nothing, on_update=nothing) where {T,V,E}
    mol = JSON.parse(json, MolGraph{T,V,E})
    isnothing(on_init) || setproperty!(mol.config, :on_init, on_init)
    isnothing(on_update) || setproperty!(mol.update, :on_update, on_update)
    return mol
end

function QueryMolGraph{T,V,E}(json::AbstractString; on_init=nothing, on_update=nothing) where {T,V,E}
    mol = JSON.parse(json, QueryMolGraph{T,V,E})
    isnothing(on_init) || setproperty!(mol.config, :on_init, on_init)
    isnothing(on_update) || setproperty!(mol.update, :on_update, on_update)
    return mol
end


# JSON auto detect

function mol_from_json(json::AbstractString
        ; on_init::Union{F,Nothing}=nothing, on_update::Union{F,Nothing}=nothing) where F
    data = JSON.parse(json)
    mol = nothing
    if haskey(data, "commonchem")
        mol = StructUtils.make(MolGraph{Int,CommonChemAtom,CommonChemBond}, data)
    elseif !haskey(data, "vproptype") || !haskey(data, "eproptype")
        error("Invalid JSON format")
    elseif data["vproptype"] == "SDFAtom" && data["eproptype"] == "SDFBond"
        mol = StructUtils.make(MolGraph{Int,SDFAtom,SDFBond}, data)
    elseif data["vproptype"] == "SMILESAtom" && data["eproptype"] == "SMILESBond"
        mol = StructUtils.make(MolGraph{Int,SMILESAtom,SMILESBond}, data)
    elseif data["vproptype"] == "QueryAtom" && data["eproptype"] == "QueryBond"
        mol = StructUtils.make(QueryMolGraph{Int,QueryAtom{Int,QueryNode},QueryBond{Int,QueryNode}}, data)
    end
    if !isnothing(mol)
        isnothing(on_init) || setproperty!(mol.config, :on_init, on_init)
        isnothing(on_update) || setproperty!(mol.update, :on_update, on_update)
        return mol
    end
    error("Invalid JSON format")
end
