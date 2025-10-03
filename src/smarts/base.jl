#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

function remap!(::Val{:smarts_lexical_succ}, gprop::SimpleMolProperty{T},
        vmap::Vector{T}, edges::Vector{Edge{T}}) where T <: Integer
    revv = Dict(v => i for (i, v) in enumerate(vmap))
    vec = [T[] for i in 1:length(revv)]
    for (k, v) in revv
        k <= length(gprop.smarts_lexical_succ) || continue
        vec[v] = [revv[s] for s in gprop.smarts_lexical_succ[k] if haskey(revv, s)]
    end
    empty!(gprop.smarts_lexical_succ)
    append!(gprop.smarts_lexical_succ, vec)
    return
end


abstract type AbstractSMARTSParser{T,V,E} end


mutable struct SMILESParser{T<:Integer,V<:AbstractAtom,E<:AbstractBond} <: AbstractSMARTSParser{T,V,E}
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

SMILESParser{T,V,E}(smiles) where {T,V,E} = SMILESParser{T,V,E}(
    smiles, 1, false, 0, 1, 1, Dict{Int,Int}(), Vector{T}[], Edge{T}[], V[], E[])



mutable struct SMARTSParser{T<:Integer,V<:QueryTree,E<:QueryTree} <: AbstractSMARTSParser{T,V,E}
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

SMARTSParser{T,V,E}(smarts) where {T,V,E} = SMARTSParser{T,V,E}(
    smarts, 1, false, 0, 1, 1, Dict{Int,Int}(), Vector{T}[], Edge{T}[], V[], E[], Vector{T}[])


function smiles_on_init!(mol::SimpleMolGraph)
    check_valence!(mol)
    stereocenter_from_smiles!(mol)
    stereobond_from_smiles!(mol)
end

function smiles_on_update!(mol::SimpleMolGraph)
    # preprocessing
    default_atom_charge!(mol)
    default_bond_order!(mol)
    kekulize!(mol)
    # recalculate bottleneck descriptors
    sssr!(mol)
    apparent_valence!(mol)
    valence!(mol)
    lone_pair!(mol)
    is_ring_aromatic!(mol)
end

function smarts_on_init!(mol::QueryMolGraph)
    stereocenter_from_smiles!(mol)
    stereobond_from_smiles!(mol)
    preprocess!(mol)
end

const SMARTSMolGraph = QueryMolGraph{Int,QueryAtom,QueryBond}


"""
    smilestomol(smiles::AbstractString) -> MolGraph{SmilesAtom,SmilesBond}

Parse SMILES string into `MolGraph` object.

The syntax of SMILES in this library follows both Daylight SMILES and OpenSMILES.

# References

1. OpenSMILES Specification http://opensmiles.org/spec/open-smiles.html
1. Daylight Tutorials https://www.daylight.com/dayhtml_tutorials/index.html
"""
function smilestomol(
        ::Type{T}, smiles::AbstractString;
        on_init=smiles_on_init!, on_update=smiles_on_update!, kwargs...) where T <: SimpleMolGraph
    V = vproptype(T) === AbstractAtom ? SMILESAtom : vproptype(T)
    E = eproptype(T) === AbstractBond ? SMILESBond : eproptype(T)
    state = SMILESParser{eltype(T),V,E}(smiles)
    fragment!(state)
    # original edge index
    gprops = MolProperty{eltype(T)}(
        smarts_lexical_succ=state.succ,
        smarts_input=string(smiles)
    )
    return T(state.edges, state.vprops, state.eprops,
        gprops=gprops, on_init=on_init, on_update=on_update; kwargs...)
end

smilestomol(smiles::AbstractString; kwargs...
    ) = smilestomol(MolGraph{Int,SMILESAtom,SMILESBond}, smiles; kwargs...)


"""
    smartstomol(smarts::AbstractString) -> QueryMolGraph{SmartsAtom,SmartsBond}

Parse SMARTS string into `QueryMolGraph` object.
"""
function smartstomol(
        ::Type{G}, smarts::AbstractString; on_init=smarts_on_init!, kwargs...) where G <: SimpleMolGraph
    T = eltype(G)
    V = vproptype(G)
    E = eproptype(G)
    state = SMARTSParser{T,V,E}(smarts)
    occursin('.', smarts) ? componentquery!(state) : fragment!(state)
    gprops = QueryMolProperty{T}(
        smarts_lexical_succ=state.succ,
        smarts_connectivity=state.connectivity,
        smarts_input=string(smarts)
    )
    return G(state.edges, state.vprops, state.eprops,
        gprops=gprops, on_init=on_init; kwargs...)
end

smartstomol(smarts::AbstractString; kwargs...
    ) = smartstomol(QueryMolGraph{Int,QueryAtom,QueryBond}, smarts; kwargs...)



# internal parser methods

function lookahead(state::AbstractSMARTSParser, pos::Int)
    # Negative pos can be used
    newpos = state.pos + pos
    if newpos > length(state.input) || newpos < 1
        return Char(0)  # \0
    else
        c = state.input[newpos]
        if isascii(c)
            return state.input[newpos]
        else
            error("smiles parser error - invalid character $(c)")
        end
    end
end

Base.read(state::AbstractSMARTSParser) = lookahead(state, 0)  # TODO: rename


function forward!(state::AbstractSMARTSParser, num::Int)
    # Negative num can be used
    state.pos += num
    if state.pos > length(state.input)
        state.done = true
    elseif state.done
        error("smiles parser error - characters in the buffer were used up")
    elseif state.pos < 1
        error("smiles parser error - no more backtracking!")
    end
end

forward!(state::AbstractSMARTSParser) = forward!(state, 1)
backtrack!(state::AbstractSMARTSParser, num::Int) = forward!(state, -num)
backtrack!(state::AbstractSMARTSParser) = backtrack!(state, 1)
