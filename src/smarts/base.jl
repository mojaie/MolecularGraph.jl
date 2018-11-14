#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    QueryMolecule,
    QueryAtom,
    QueryBond,
    AbstractSmartsParserState,
    SmartsParserState,
    SmilesParserState,
    lookahead,
    read,
    forward!,
    backtrack!,
    parse,
    SMILES,
    SMARTS


import Base: read, parse


struct QueryAtom <: AbstractNode
    query::Pair
end


struct QueryBond <: AbstractEdge
    u::Int
    v::Int
    query::Pair
end

QueryBond(b::QueryBond, u, v) = QueryBond(u, v, b.query)


struct QueryMolecule <: AbstractMutableMolecule
    graph::MutableUDGraph{QueryAtom,QueryBond}
    connectivity::Array{Array{Int}}
    attribute::Dict

    function QueryMolecule()
        new(MutableUDGraph{QueryAtom,QueryBond}(), [], Dict())
    end
end


abstract type AbstractSmartsParserState end


mutable struct SmilesParserState <: AbstractSmartsParserState
    pos::Int
    done::Bool
    input::String
    ringlabel::Dict
    mol::MutableMolecule

    function SmilesParserState(smiles)
        new(1, false, smiles, Dict(), MutableMolecule())
    end
end


mutable struct SmartsParserState <: AbstractSmartsParserState
    pos::Int
    done::Bool
    input::String
    ringlabel::Dict
    mol::QueryMolecule

    function SmartsParserState(smarts)
        new(1, false, smarts, Dict(), QueryMolecule())
    end
end


SMILES = MutableMolecule
SMARTS = QueryMolecule


function parse(::Type{T}, str::AbstractString) where T <: SMILES
    state = SmilesParserState(str)
    connectedquery!(state)
    state.mol.attribute[:sourcetype] = :smiles
    return state.mol
end


function parse(::Type{T}, str::AbstractString) where T <: SMARTS
    state = SmartsParserState(str)
    componentquery!(state)
    state.mol.attribute[:sourcetype] = :smarts
    return state.mol
end


function lookahead(state::AbstractSmartsParserState, pos::Int)
    # Negative pos can be used
    newpos = state.pos + pos
    if newpos > length(state.input) || newpos < 1
        return Char(0)
    else
        c = state.input[newpos]
        if isascii(c)
            return state.input[newpos]
        else
            throw(ParserError("invalid charactor $(c)"))
        end
    end
end

read(state) = lookahead(state, 0)


function forward!(state::AbstractSmartsParserState, num::Int)
    # Negative num can be used
    state.pos += num
    if state.pos > length(state.input)
        state.done = true
    elseif state.done
        throw(ParserError("Charactors in the buffer were used up"))
    elseif state.pos < 1
        throw(ParserError("No more backtracking!"))
    end
end

forward!(state) = forward!(state, 1)
backtrack!(state, num) = forward!(state, -num)
backtrack!(state) = backtrack!(state, 1)
