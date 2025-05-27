#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    SMILESParser, SMARTSParser,
    smilestomol, smartstomol


struct SMARTSBondIndex{T} <: AbstractDict{Edge{T},T}
    mapping::Dict{Edge{T},T}
end

SMARTSBondIndex{T}(data::Vector=[]) where T = SMARTSBondIndex{T}(
    Dict{Edge{T},T}(u_edge(T, e[1], e[2]) => i for (e, i) in data))

Base.iterate(bi::SMARTSBondIndex) = iterate(bi.mapping)
Base.iterate(bi::SMARTSBondIndex, i) = iterate(bi.mapping, i)
Base.length(bi::SMARTSBondIndex) = length(bi.mapping)
Base.get(bi::SMARTSBondIndex, k, v) = get(bi.mapping, k, v)
Base.setindex!(bi::SMARTSBondIndex, v, k) = setindex!(bi.mapping, v, k)
to_dict(::Val{:default}, bi::SMARTSBondIndex) = [[[src(e), dst(e)], i] for (e, i) in bi.mapping]

function remap(bi::SMARTSBondIndex{T}, vmap::Dict) where T  # vmap[old] -> new
    newmap = Dict{Edge{T},T}()
    for (e, i) in bi.mapping
        (haskey(vmap, src(e)) && haskey(vmap, dst(e))) || continue
        newmap[u_edge(T, vmap[src(e)], vmap[dst(e)])] = i
    end
    return SMARTSBondIndex{T}(newmap)
end


mutable struct SMILESParser{T,V,E}
    input::String
    pos::Int
    done::Bool
    node::Int # No. of current node
    branch::Int # No. of node at the current branch root
    root::Int # No. of node at the current tree root
    ringlabel::Dict{Int,Int}  # label => original bond index
    edges::Vector{Edge{T}}
    vprops::Vector{V}
    eprops::Vector{E}
end

SMILESParser{T}(smiles
    ) where T <: SimpleMolGraph = SMILESParser{eltype(T),vproptype(T),eproptype(T)}(
        smiles, 1, false, 0, 1, 1, Dict(), Edge{eltype(T)}[], vproptype(T)[], eproptype(T)[])



mutable struct SMARTSParser{T,V,E}
    input::String
    pos::Int
    done::Bool
    node::Int # No. of current node
    branch::Int # No. of node at the current branch root
    root::Int # No. of node at the current tree root
    ringlabel::Dict{Int,Int}  # label => original bond index
    edges::Vector{Edge{T}}
    vprops::Vector{V}
    eprops::Vector{E}
    connectivity::Vector{Vector{T}}
end

SMARTSParser{T}(smarts
    ) where T <: SimpleMolGraph = SMARTSParser{eltype(T),vproptype(T),eproptype(T)}(
        smarts, 1, false, 0, 1, 1, Dict(), Edge{eltype(T)}[], vproptype(T)[], eproptype(T)[], [])


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
        :original_bond_index => SMARTSBondIndex{eltype(T)}(
            Dict(e => i for (i, e) in enumerate(state.edges))),
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
        :original_bond_index => SMARTSBondIndex{eltype(T)}(
            Dict(e => i for (i, e) in enumerate(state.edges))),
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
