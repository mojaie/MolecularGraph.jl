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


"""
    SMILESParser{T<:Integer,V<:AbstractAtom,E<:AbstractBond}

Parser state of a SMILES string used in [`smilestomol`](@ref).
"""
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
    bracket::Vector{T}  # nodes of bracket atoms (e.g. [CH2])
end

SMILESParser{T,V,E}(smiles) where {T,V,E} = SMILESParser{T,V,E}(
    smiles, 1, false, 0, 1, 1, Dict{Int,Int}(), Vector{T}[], Edge{T}[], V[], E[], T[])



"""
    SMARTSParser{T<:Integer,V<:QueryTree,E<:QueryTree}

Parser state of a SMARTS string used in [`smartstomol`](@ref).
"""
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


"""
    smiles_on_init!(mol::SimpleMolGraph) -> Nothing

Default `on_init` callback of [`smilestomol`](@ref) (valence check and stereochemistry).
"""
function smiles_on_init!(mol::SimpleMolGraph)
    check_valence!(mol)
    stereocenter_from_smiles!(mol)
    stereobond_from_smiles!(mol)
end

"""
    smiles_on_update!(mol::SimpleMolGraph) -> Nothing

Default `on_update` callback of [`smilestomol`](@ref) (standardization of charges and bond
orders, kekulization and recalculation of cached descriptors).
"""
function smiles_on_update!(mol::SimpleMolGraph)
    # preprocessing
    default_atom_charge!(mol)
    default_bond_order!(mol)
    kekulize!(mol)  # also sets apparent_valence
    # recalculate bottleneck descriptors
    sssr!(mol)
    valence!(mol)
    lone_pair!(mol)
    is_ring_aromatic!(mol)
end

"""
    smarts_on_init!(mol::QueryMolGraph) -> Nothing

Default `on_init` callback of [`smartstomol`](@ref) (stereochemistry and query preprocessing).
"""
function smarts_on_init!(mol::QueryMolGraph)
    stereocenter_from_smiles!(mol)
    stereobond_from_smiles!(mol)
    preprocess!(mol)
end

"""
    SMARTSMolGraph

Alias of `QueryMolGraph{Int,QueryAtom{Int,QueryNode},QueryBond{Int,QueryNode}}`,
the molecule type of [`smartstomol`](@ref).
"""
const SMARTSMolGraph = QueryMolGraph{Int,QueryAtom{Int,QueryNode},QueryBond{Int,QueryNode}}


"""
    smilestomol(smiles::AbstractString; kwargs...) -> MolGraph{Int,SMILESAtom,SMILESBond}
    smilestomol(::Type{T}, smiles::AbstractString; on_init=smiles_on_init!, on_update=smiles_on_update!) -> T

Parse SMILES string into `MolGraph` object.

The syntax of SMILES in this library follows both Daylight SMILES and OpenSMILES.

# References

1. OpenSMILES Specification http://opensmiles.org/spec/open-smiles.html
1. Daylight Tutorials https://www.daylight.com/dayhtml_tutorials/index.html
"""
function smilestomol(
        ::Type{T}, smiles::AbstractString;
        on_init::F1=smiles_on_init!, on_update::F2=smiles_on_update!, kwargs...
        ) where {T<:SimpleMolGraph,F1,F2}  # F1, F2: force specialization on callbacks
    V = vproptype(T) === AbstractAtom ? SMILESAtom : vproptype(T)
    E = eproptype(T) === AbstractBond ? SMILESBond : eproptype(T)
    state = SMILESParser{eltype(T),V,E}(smiles)
    fragment!(state)
    bracket_valence!(state)
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


# Bracket atoms have no implicit hydrogens (e.g. [CH2] is a carbene, [CH3] is a methyl radical).
# If the total bond order (including hydrogens) is lower than the typical valence, it is set
# to the atom as the specified valence. The total bond order of aromatic atoms is the same
# as `apparent_valence` before kekulization (aromatic bonds as single bonds plus a pi bond,
# see `aromatic_pi_bond`), e.g. [c] is an aryl radical with valence 3.
function bracket_valence!(state::SMILESParser{T,V,E}) where {T,V,E}
    isempty(state.bracket) && return
    # hydrogens are graph nodes (hydrogen count 0 for bracket atoms)
    hcounts = fill(-1, length(state.vprops))
    hcounts[state.bracket] .= 0
    arom = BitVector(has_isaromatic(V) && isaromatic(a) for a in state.vprops)
    vals = hydrogen_count_valence(
        state.edges, Int[bond_order(b) for b in state.eprops], arom,
        Symbol[atom_symbol(a) for a in state.vprops], Int[atom_charge(a) for a in state.vprops],
        hcounts)
    for i in state.bracket
        vals[i] == -1 && continue
        a = state.vprops[i]
        state.vprops[i] = arom[i] ? V(; symbol=atom_symbol(a), charge=atom_charge(a),
                isotope=isotope(a), isaromatic=true, stereo=a.stereo, valence=vals[i]) :
            V(; symbol=atom_symbol(a), charge=atom_charge(a), isotope=isotope(a),
                stereo=a.stereo, valence=vals[i])
    end
    return
end


"""
    smartstomol(smarts::AbstractString; kwargs...) -> SMARTSMolGraph
    smartstomol(::Type{G}, smarts::AbstractString; on_init=smarts_on_init!, kwargs...) -> G

Parse SMARTS string into `QueryMolGraph` object.
"""
function smartstomol(
        ::Type{G}, smarts::AbstractString; on_init::F=smarts_on_init!, kwargs...
        ) where {G<:SimpleMolGraph,F}  # F: force specialization on callbacks
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
    ) = smartstomol(SMARTSMolGraph, smarts; kwargs...)



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
