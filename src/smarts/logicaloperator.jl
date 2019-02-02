#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    lglowand!,
    lgor!,
    lghighand!,
    lgnot!

""" The argument `func` is a parser function which has a parser state as an
argument, process tokens found in the given text, and returns nothing if no
valid tokens were found.
"""

function lglowand!(state::AnySmarts, func)
    """ LogicalLowAnd <- Or (';' Or)*
    """
    cond = []
    c = lgor!(state, func)
    while c !== nothing
        push!(cond, c)
        if read(state) == ';'
            forward!(state)
            c = lgor!(state, func)
        else
            break
        end
    end
    if isempty(cond)
        return
    elseif length(cond) == 1
        return cond[1]
    else
        return :and => Tuple(cond)
    end
end


function lgor!(state::AnySmarts, func)
    """ Or <- And (',' And)*
    """
    cond = []
    c = lghighand!(state, func)
    while c !== nothing
        push!(cond, c)
        if read(state) == ','
            forward!(state)
            c = lghighand!(state, func)
        else
            break
        end
    end
    if isempty(cond)
        return
    elseif length(cond) == 1
        return cond[1]
    else
        return :or => Tuple(cond)
    end
end


function lghighand!(state::SmartsParser, func)
    """ And <- Not ('&'? Not)*
    """
    cond = []
    c = lgnot!(state, func)
    while c !== nothing
        push!(cond, c)
        if read(state) == '&'
            forward!(state)
        end
        c = lgnot!(state, func)
    end
    if isempty(cond)
        return
    elseif length(cond) == 1
        return cond[1]
    else
        return :and => Tuple(cond)
    end
end


function lgnot!(state::SmartsParser, func)
    """ Not <- '!'? Element
    """
    if read(state) == '!'
        forward!(state)
        cond = func(state)
        if cond === nothing
            # TODO: need backtracking in the case like [!=!c]
            backtrack!(state)
            return
        else
            return :not => cond
        end
    else
        return func(state)
    end
end
