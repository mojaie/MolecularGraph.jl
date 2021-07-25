#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    SmilesParser, SmartsParser,
    smilestomol, smartstomol,
    associate_operations


mutable struct SmartsParserState{N<:AbstractNode,E<:UndirectedEdge}
    input::String
    allow_disconnected::Bool

    pos::Int
    done::Bool
    node::Int # No. of current node
    branch::Int # No. of node at the current branch root
    root::Int # No. of node at the current tree root
    ringlabel::Dict{Int,Tuple{Int,Union{E,Symbol,Nothing}}} # TODO: strange union

    edges::Vector{Tuple{Int,Int}}
    nodeattrs::Vector{N}
    edgeattrs::Vector{E}
    connectivity::Vector{Vector{Int}}

    function SmartsParserState{N,E}(str, disconn
            ) where {N<:AbstractNode,E<:UndirectedEdge}
        new(str, disconn, 1, false, 0, 1, 1, Dict(), [], [], [], [])
    end
end

SmilesParser = SmartsParserState{SmilesAtom,SmilesBond}
SmartsParser = SmartsParserState{SmartsAtom,SmartsBond}


function Base.parse(::Type{SMILES}, str::AbstractString)
    state = SmilesParser(str, true)
    fragment!(state)
    return graphmol(state.edges, state.nodeattrs, state.edgeattrs)
end

function Base.parse(::Type{SMARTS}, str::AbstractString)
    if occursin('.', str)
        state = SmartsParser(str, true)
        componentquery!(state)
    else
        state = SmartsParser(str, false)
        fragment!(state)
    end
    qmol = querymol(state.edges, state.nodeattrs, state.edgeattrs, state.connectivity)
    resolvedefaultbond!(qmol)
    return qmol
end


"""
    smilestomol(smiles::AbstractString) -> GraphMol{SmilesAtom,SmilesBond}

Parse SMILES string into `GraphMol` object.

The syntax of SMILES in this library follows both Daylight SMILES and OpenSMILES.

# References

1. OpenSMILES Specification http://opensmiles.org/spec/open-smiles.html
1. Daylight Tutorials https://www.daylight.com/dayhtml_tutorials/index.html
"""
function smilestomol(smiles::AbstractString)
    mol = parse(SMILES, smiles)
    kekulize!(mol)
    setdiastereo!(mol)
    return mol
end


"""
    smartstomol(smarts::AbstractString) -> QueryMol{SmartsAtom,SmartsBond}

Parse SMARTS string into `QueryMol` object.
"""
smartstomol(smarts::AbstractString) = parse(SMARTS, smarts)



function lookahead(state::SmartsParserState, pos::Int)
    # Negative pos can be used
    newpos = state.pos + pos
    if newpos > length(state.input) || newpos < 1
        return Char(0)
    else
        c = state.input[newpos]
        if isascii(c)
            return state.input[newpos]
        else
            throw(ErrorException("invalid charactor $(c)"))
        end
    end
end

Base.read(state) = lookahead(state, 0)  # TODO: rename


function forward!(state::SmartsParserState, num::Int)
    # Negative num can be used
    state.pos += num
    if state.pos > length(state.input)
        state.done = true
    elseif state.done
        throw(ErrorException("charactors in the buffer were used up"))
    elseif state.pos < 1
        throw(ErrorException("no more backtracking!"))
    end
end

forward!(state) = forward!(state, 1)
backtrack!(state, num) = forward!(state, -num)
backtrack!(state) = backtrack!(state, 1)


"""
    associate_operations(fml::Pair) -> Pair

Return juxtaposed formulae (ex. :and => (A, :and => (B, C)) -> :and => (A, B, C)).
"""
function associate_operations(fml::Pair)
    fml.first in (:and, :or) || return fml
    associated = Set()
    for elem in fml.second
        if elem.first === fml.first
            union!(associated, associate_operations(elem).second)
        else
            push!(associated, elem)
        end
    end
    return fml.first => Tuple(associated)
end


"""
    resolvedefaultbond(qmol::QueryMol) -> Nothing

Resolve default SMARTS bonds.
"""
function resolvedefaultbond!(qmol::QueryMol)
    newedges = QueryBond[]
    for e in 1:edgecount(qmol)
        (u, v) = getedge(qmol, e)
        q = edgeattr(qmol, e).query
        if q == (:defaultbond => true)
            uq = nodeattr(qmol, u).query
            isuqalip = uq == (:isaromatic => false) || (uq.first === :and && (:isaromatic => false) in uq.second)
            isuqarom = uq == (:isaromatic => true) || (uq.first === :and && (:isaromatic => true) in uq.second)
            vq = nodeattr(qmol, v).query
            isvqalip = vq == (:isaromatic => false) || (uq.first === :and && (:isaromatic => false) in uq.second)
            isvqarom = vq == (:isaromatic => true) || (vq.first === :and && (:isaromatic => true) in vq.second)
            if isuqarom && isvqarom
                push!(newedges, SmartsBond(:isaromaticbond => true))
            elseif (isuqalip || isuqarom) && (isvqalip || isvqarom)
                push!(newedges, SmartsBond(:bondorder => 1))
            else
                push!(newedges, SmartsBond(:or => (:bondorder => 1, :isaromaticbond => true)))
            end
        else
            push!(newedges, edgeattr(qmol, e))
        end
    end
    empty!(qmol.edgeattrs)
    append!(qmol.edgeattrs, newedges)
end