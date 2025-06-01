#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

struct SMARTSLexicalSuccessors{T} <: AbstractVector{AbstractVector{T}}
    vector::Vector{Vector{T}}
end

Base.size(arr::SMARTSLexicalSuccessors) = size(arr.vector)
Base.getindex(arr::SMARTSLexicalSuccessors, i...) = getindex(arr.vector, i...)
to_dict(::Val{:default}, key::Symbol, arr::SMARTSLexicalSuccessors) = Dict{String,Any}(
    "key" => string(key),
    "type" => "SMARTSLexicalSuccessors",
    "data" => arr.vector
)
PROPERTY_TYPE_REGISTRY["SMARTSLexicalSuccessors"] = (T, data) -> SMARTSLexicalSuccessors{eltype(T)}(data)

function remap(arr::SMARTSLexicalSuccessors{T}, vmap::Dict) where T
    vec = [T[] for i in 1:length(vmap)]
    for (k, v) in vmap
        k in arr.vector || continue
        vec[v] = [vmap[s] for s in arr.vector[k] if haskey(vmap, s)]
    end
    return SMARTSLexicalSuccessors{T}(vec)
end


mutable struct SMILESParser{T,V,E}
    input::String
    pos::Int
    done::Bool
    node::Int # No. of current node
    branch::Int # No. of node at the current branch root
    root::Int # No. of node at the current tree root
    ringlabel::Dict{Int,Int}  # ring label => No. of edge to connect
    succ::Vector{Vector{T}}  # node => successors in lexical order (for stereochem)
    edges::Vector{Edge{T}}
    vprops::Vector{V}
    eprops::Vector{E}
end

SMILESParser{T}(smiles
    ) where T <: SimpleMolGraph = SMILESParser{eltype(T),vproptype(T),eproptype(T)}(
        smiles, 1, false, 0, 1, 1, Dict(), [], Edge{eltype(T)}[], vproptype(T)[], eproptype(T)[])



mutable struct SMARTSParser{T,V,E}
    input::String
    pos::Int
    done::Bool
    node::Int # No. of current node
    branch::Int # No. of node at the current branch root
    root::Int # No. of node at the current tree root
    ringlabel::Dict{Int,Int}  # label => original bond index
    succ::Vector{Vector{T}}  # node => successors in lexical order (for stereochem)
    edges::Vector{Edge{T}}
    vprops::Vector{V}
    eprops::Vector{E}
    connectivity::Vector{Vector{T}}
end

SMARTSParser{T}(smarts
    ) where T <: SimpleMolGraph = SMARTSParser{eltype(T),vproptype(T),eproptype(T)}(
        smarts, 1, false, 0, 1, 1, Dict(), [], Edge{eltype(T)}[], vproptype(T)[], eproptype(T)[], [])


function smiles_on_init!(mol)
    update_edge_rank!(mol)
    stereocenter_from_smiles!(mol)
    stereobond_from_smiles!(mol)
    set_state!(mol, :initialized, true)
end

function smiles_on_update!(mol)
    update_edge_rank!(mol)
    clear_caches!(mol)
    set_state!(mol, :has_updates, false)
    # preprocessing
    kekulize!(mol)
    # recalculate bottleneck descriptors
    sssr!(mol)
    lone_pair!(mol)
    apparent_valence!(mol)
    valence!(mol)
    is_ring_aromatic!(mol)
end


"""
    smilestomol(smiles::AbstractString) -> GraphMol{SmilesAtom,SmilesBond}

Parse SMILES string into `GraphMol` object.

The syntax of SMILES in this library follows both Daylight SMILES and OpenSMILES.

# References

1. OpenSMILES Specification http://opensmiles.org/spec/open-smiles.html
1. Daylight Tutorials https://www.daylight.com/dayhtml_tutorials/index.html
"""
function smilestomol(::Type{T}, smiles::AbstractString;
        config=Dict{Symbol,Any}(), kwargs...) where T <: AbstractMolGraph
    state = SMILESParser{T}(smiles)
    fragment!(state)
    # original edge index
    gprops = Dict(
        :lexical_successors => SMARTSLexicalSuccessors{eltype(T)}(state.succ),
        :metadata => Metadata()
    )
    default_config = Dict{Symbol,Any}(:updater => smiles_on_update!, :on_init => smiles_on_init!)
    merge!(default_config, config)
    return T(state.edges, state.vprops, state.eprops,
        gprop_map=gprops, config_map=default_config; kwargs...)
end

smilestomol(smiles::AbstractString; kwargs...
    ) = smilestomol(MolGraph{Int,SMILESAtom,SMILESBond}, smiles; kwargs...)


"""
    smartstomol(smarts::AbstractString) -> QueryMol{SmartsAtom,SmartsBond}

Parse SMARTS string into `QueryMol` object.
"""
function smartstomol(::Type{T}, smarts::AbstractString;
        gprop_map=Dict{Symbol,Any}(), kwargs...) where T <: AbstractMolGraph
    state = SMARTSParser{T}(smarts)
    occursin('.', smarts) ? componentquery!(state) : fragment!(state)
    default_gprop = Dict{Symbol,Any}(
        :connectivity => state.connectivity,
        :lexical_successors => SMARTSLexicalSuccessors{eltype(T)}(state.succ),
        :metadata => Metadata()
    )
    merge!(default_gprop, gprop_map)
    mol = T(
        state.edges, state.vprops, state.eprops, gprop_map=default_gprop; kwargs...)
    if vproptype(mol) <: QueryTree  # vproptype(mol) can be QueryTruthTable for testing
        specialize_nonaromatic!(mol)
        remove_hydrogens!(mol)
        for i in vertices(mol)
            set_prop!(mol, i, QueryTree(optimize_query(get_prop(mol, i, :tree))))
        end
    end
    return mol
end

smartstomol(smarts::AbstractString
    ) = smartstomol(MolGraph{Int,QueryTree,QueryTree}, smarts)



# internal parser methods

function lookahead(state::Union{SMILESParser,SMARTSParser}, pos::Int)
    # Negative pos can be used
    newpos = state.pos + pos
    if newpos > length(state.input) || newpos < 1
        return Char(0)
    else
        c = state.input[newpos]
        if isascii(c)
            return state.input[newpos]
        else
            error("invalid charactor $(c)")
        end
    end
end

Base.read(state) = lookahead(state, 0)  # TODO: rename


function forward!(state::Union{SMILESParser,SMARTSParser}, num::Int)
    # Negative num can be used
    state.pos += num
    if state.pos > length(state.input)
        state.done = true
    elseif state.done
        error("charactors in the buffer were used up")
    elseif state.pos < 1
        error("no more backtracking!")
    end
end

forward!(state) = forward!(state, 1)
backtrack!(state, num) = forward!(state, -num)
backtrack!(state) = backtrack!(state, 1)
