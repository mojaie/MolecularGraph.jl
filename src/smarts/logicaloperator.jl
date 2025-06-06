#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#


"""
    lglowand!(state::SmartsParserState) -> Union{Pair,Nothing}

LogicalLowAnd <- Or (';' Or)*

The argument `func` is a parser function which has a parser state as an argument,
process tokens found in the given text, and returns nothing if no valid tokens were found.
"""
function lglowand!(state::SMARTSParser, func)
    qs = []
    q = lgor!(state, func)
    q isa EndToken && return EndToken()
    while !isa(q, EndToken)
        push!(qs, q)
        if readtoken(state) == ';'
            forward!(state)
            q = lgor!(state, func)
            continue
        end
        break
    end
    isempty(qs) && error("(lglowand!) invalid AND(;) operation")
    length(qs) == 1 && return qs[1]
    return QueryOperator(:and, qs)
end


"""
    lgor!(state::SmartsParserState) -> Union{Pair,Nothing}

Or <- And (',' And)*

The argument `func` is a parser function which has a parser state as an argument,
process tokens found in the given text, and returns nothing if no valid tokens were found.
"""
function lgor!(state::SMARTSParser, func)
    qs = []
    q = lghighand!(state, func)
    q isa EndToken && return EndToken()
    while !isa(q, EndToken)
        push!(qs, q)
        if readtoken(state) == ','
            forward!(state)
            q = lghighand!(state, func)
            continue
        end
        break
    end
    isempty(qs) && error("(lgor!) invalid OR(,) operation")
    length(qs) == 1 && return qs[1]
    return QueryOperator(:or, qs)
end


"""
    lghighand!(state::SmartsParserState) -> Union{Pair,Nothing}

And <- Not ('&'? Not)*

The argument `func` is a parser function which has a parser state as an argument,
process tokens found in the given text, and returns nothing if no valid tokens were found.
"""
function lghighand!(state::T, func) where T <: AbstractSMARTSParser
    qs = QueryComponent[]
    q = gethighand(state, func)
    q isa EndToken && return EndToken()
    while !isa(q, EndToken)
        push!(qs, q)
        readtoken(state) == '&' && forward!(state)
        q = gethighand(state, func)
    end
    isempty(qs) && error("(lghighand!) invalid AND(&) operation")
    length(qs) == 1 && return qs[1]
    return QueryOperator(:and, qs)
end

gethighand(state::SMILESParser, func) = func(state)
gethighand(state::SMARTSParser, func) = lgnot!(state, func)


"""
    lgnot!(state::SmartsParserState) -> Union{Pair,Nothing}

Not <- '!'? Element

The argument `func` is a parser function which has a parser state as an argument,
process tokens found in the given text, and returns nothing if no valid tokens were found.
"""
function lgnot!(state::SMARTSParser, func)
    if readtoken(state) == '!'
        forward!(state)
        q = func(state)
        q isa EndToken && error("(lgnot!) invalid NOT(!) operation")
        return QueryOperator(:not, [q])
    end
    return func(state)  # can return EndToken if the parser get stop token
end
