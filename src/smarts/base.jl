#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    SmartsParser,
    parse,
    smilestomol,
    lookahead,
    read,
    forward!,
    backtrack!


import Base: read, parse


mutable struct SmartsParser{T<:Union{VectorMol,QueryMol}}
    input::String
    allow_disconnected::Bool

    pos::Int
    done::Bool
    node::Int # No. of current node
    branch::Int # No. of node at the current branch root
    root::Int # No. of node at the current tree root
    ringlabel::Dict

    mol::T

    function SmartsParser{T}(str, disconn) where {T<:Union{VectorMol,QueryMol}}
        new(str, disconn, 1, false, 0, 1, 1, Dict(), T())
    end
end


function parse(::Type{SMILES}, str::AbstractString)
    state = SmartsParser{SMILES}(str, true)
    fragment!(state)
    return state.mol
end

function parse(::Type{SMARTS}, str::AbstractString)
    if occursin('.', str)
        state = SmartsParser{SMARTS}(str, true)
        componentquery!(state)
    else
        state = SmartsParser{SMARTS}(str, false)
        fragment!(state)
    end
    return state.mol
end


function smilestomol(smiles::AbstractString)
    return parse(SMILES, smiles)
end


function lookahead(state::SmartsParser, pos::Int)
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

read(state) = lookahead(state, 0)


function forward!(state::SmartsParser, num::Int)
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
