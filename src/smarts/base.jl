#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    AbstractSmartsParser,
    SmartsParserState,
    SmilesParserState,
    parse,
    smilestomol,
    lookahead,
    read,
    forward!,
    backtrack!


import Base: read, parse


abstract type AbstractSmartsParser end


mutable struct SmilesParserState <: AbstractSmartsParser
    input::String
    pos::Int
    done::Bool
    node::Int # No. of current node
    branch::Int # No. of node at the current branch root
    root::Int # No. of node at the current tree root
    ringlabel::Dict
    mol::SMILES

    function SmilesParserState(smiles)
        new(smiles, 1, false, 0, 1, 1, Dict(), SMILES())
    end
end


mutable struct SmartsParserState <: AbstractSmartsParser
    input::String
    pos::Int
    done::Bool
    node::Int # No. of current node
    branch::Int # No. of node at the current branch root
    root::Int # No. of node at the current tree root
    ringlabel::Dict
    mol::SMARTS

    function SmartsParserState(smarts)
        new(smarts, 1, false, 0, 1, 1, Dict(), SMARTS())
    end
end


function parse(::Type{T}, str::AbstractString) where T <: SMILES
    state = SmilesParserState(str)
    fragment!(state)
    return state.mol
end


function parse(::Type{T}, str::AbstractString) where T <: SMARTS
    state = SmartsParserState(str)
    componentquery!(state)
    return state.mol
end


function smilestomol(smiles::AbstractString)
    mol = parse(SMILES, smiles)
    vmol = vectormol(mol)
    default_annotation!(vmol)
    return vmol
end


function lookahead(state::AbstractSmartsParser, pos::Int)
    # Negative pos can be used
    newpos = state.pos + pos
    if newpos > length(state.input) || newpos < 1
        return Char(0)
    else
        c = state.input[newpos]
        if isascii(c)
            return state.input[newpos]
        else
            throw(MolParseError("invalid charactor $(c)"))
        end
    end
end

read(state) = lookahead(state, 0)


function forward!(state::AbstractSmartsParser, num::Int)
    # Negative num can be used
    state.pos += num
    if state.pos > length(state.input)
        state.done = true
    elseif state.done
        throw(MolParseError("Charactors in the buffer were used up"))
    elseif state.pos < 1
        throw(MolParseError("No more backtracking!"))
    end
end

forward!(state) = forward!(state, 1)
backtrack!(state, num) = forward!(state, -num)
backtrack!(state) = backtrack!(state, 1)
