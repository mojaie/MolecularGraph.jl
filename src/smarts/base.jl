#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    SmartsAtom, SmartsBond, SMARTS,
    SmilesParser, SmartsParser,
    smilestomol, smartstomol,
    associate_operations


struct SmartsAtom <: QueryAtom
    query::QueryFormula
end


struct SmartsBond <: QueryBond
    query::QueryFormula
end

SmartsBond() = SmartsBond(QueryFormula(:any, true))


SMARTS = QueryMol{SmartsAtom,SmartsBond}



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
    resolvedefaultbond(qmol::QueryMol) -> Nothing

Resolve default SMARTS bonds.
"""
function resolvedefaultbond!(qmol::QueryMol)
    for e in 1:edgecount(qmol)
        q = edgeattr(qmol, e).query
        q == QueryFormula(:defaultbond, true) || continue
        (u, v) = getedge(qmol, e)
        uq = nodeattr(qmol, u).query
        uarom = findformula(uq, :isaromatic)
        vq = nodeattr(qmol, v).query
        varom = findformula(vq, :isaromatic)
        arombond = QueryFormula(:isaromaticbond, true)
        single = QueryFormula(:bondorder, 1)
        if uarom === true && varom === true
            setedgeattr!(qmol, e, SmartsBond(arombond))
        elseif uarom === false || varom === false
            setedgeattr!(qmol, e, SmartsBond(single))
        else
            setedgeattr!(qmol, e, SmartsBond(QueryFormula(:or, Set([single, arombond]))))
        end
    end
end
