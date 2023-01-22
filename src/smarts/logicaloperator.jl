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
    q === nothing && return
    while q !== nothing
        push!(qs, q)
        if read(state) == ';'
            forward!(state)
            q = lgor!(state, func)
            continue
        end
        break
    end
    @assert !isempty(qs) "(lglowand!) invalid AND(;) operation"
    length(qs) == 1 && return qs[1]
    return and_query(qs)
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
    q === nothing && return
    while q !== nothing
        push!(qs, q)
        if read(state) == ','
            forward!(state)
            q = lghighand!(state, func)
            continue
        end
        break
    end
    @assert !isempty(qs) "(lgor!) invalid OR(,) operation"
    length(qs) == 1 && return qs[1]
    return or_query(qs)
end


"""
    lghighand!(state::SmartsParserState) -> Union{Pair,Nothing}

And <- Not ('&'? Not)*

The argument `func` is a parser function which has a parser state as an argument,
process tokens found in the given text, and returns nothing if no valid tokens were found.
"""
function lghighand!(state::Union{SMILESParser,SMARTSParser}, func)
    qs = []
    q = isa(state, SMILESParser) ? func(state) : lgnot!(state, func)
    q === nothing && return
    while q !== nothing
        push!(qs, q)
        read(state) == '&' && forward!(state)
        q = isa(state, SMILESParser) ? func(state) : lgnot!(state, func)
    end
    @assert !isempty(qs) "(lghighand!) invalid AND(&) operation"
    length(qs) == 1 && return qs[1]
    return and_query(qs)
end


"""
    lgnot!(state::SmartsParserState) -> Union{Pair,Nothing}

Not <- '!'? Element

The argument `func` is a parser function which has a parser state as an argument,
process tokens found in the given text, and returns nothing if no valid tokens were found.
"""
function lgnot!(state::SMARTSParser, func)
    if read(state) == '!'
        forward!(state)
        q = func(state)
        @assert q !== nothing "(lgnot!) invalid NOT(!) operation"
        return not_query(q)
    end
    return func(state)  # can be Nothing if the parser get stop token
end
