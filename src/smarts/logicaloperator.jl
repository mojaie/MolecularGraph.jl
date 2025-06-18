#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#


function lgoperator!(
        state::AbstractSMARTSParser, qtree::QueryTree,
        token::Char, downstream::F, op::QueryNode) where F
    vs = Int[]
    v = downstream(state, qtree)
    v == 0 && return 0  # not found
    while v != 0
        push!(vs, v)
        if read(state) == token
            forward!(state)
            v = downstream(state, qtree)
        elseif token == '&'
            v = downstream(state, qtree)
        else
            break
        end
    end
    isempty(vs) && error("($(string(downstream))) invalid operator '$(string(token))'")
    length(vs) == 1 && return vs[1]
    node = add_qnode!(qtree, op)
    for i in vs
        add_qedge!(qtree, node, i)
    end
    return node
end


"""
    lglowand!(state::SmartsParser, qtree::T) where {T<:QueryTree} -> Nothing

LogicalLowAnd <- Or (';' Or)*

The argument `func` is a parser function which has a parser state as an argument,
process tokens found in the given text, and returns nothing if no valid tokens were found.
"""
lglowand!(state::SMARTSParser, qtree::QueryTree
) = lgoperator!(state, qtree, ';', lgor!, qand())


"""
    lgor!(state::SmartsParserState) -> Union{Pair,Nothing}

Or <- And (',' And)*

The argument `func` is a parser function which has a parser state as an argument,
process tokens found in the given text, and returns nothing if no valid tokens were found.
"""
lgor!(state::SMARTSParser, qtree::QueryTree
) = lgoperator!(state, qtree, ',', lghighand!, qor())


"""
    lghighand!(state::SmartsParserState) -> Union{Pair,Nothing}

And <- Not ('&'? Not)*

The argument `func` is a parser function which has a parser state as an argument,
process tokens found in the given text, and returns nothing if no valid tokens were found.
"""
lghighand!(state::SMARTSParser, qtree::QueryTree
    ) = lgoperator!(state, qtree, '&', lgnot!, qand())

function lghighand!(state::SMILESParser{T,V,E}, downstream::F) where {T,V,E,F}
    vs = NamedTuple[]
    hcnt = 0
    v = downstream(state)
    isnothing(v) && return  # not found
    while !isnothing(v)
        if haskey(v, :total_hydrogens)
            hcnt = v[:total_hydrogens]
        else
            push!(vs, v)
        end
        v = downstream(state)
    end
    isempty(vs) && hcnt == 0 && error("($(string(downstream))) invalid operator '&'")
    # explicit hydrogen (e.g. [CH3]) -> hydrogen nodes
    merged = isempty(vs) ? NamedTuple() : merge(vs...)
    if !haskey(merged, :symbol)  # special case: [H2]
        merged = merge(merged, (symbol=:H,))
        hcnt -= 1
    end
    return [V(;merged...), (V(;symbol=:H) for _ in 1:hcnt)...]
end

"""
    lgnot!(state::SmartsParserState) -> Union{Pair,Nothing}

Not <- '!'? Element

The argument `func` is a parser function which has a parser state as an argument,
process tokens found in the given text, and returns nothing if no valid tokens were found.
"""
function lgnot!(state::SMARTSParser, qtree::QueryTree, downstream::F) where F
    if read(state) == '!'
        forward!(state)
        v = downstream(state, qtree)
        v == 0 && error("(lgnot!) invalid NOT(!) operation")
        node = add_qnode!(qtree, qnot())
        add_qedge!(qtree, node, v)
        return node
    else
        v = downstream(state, qtree)
        v == 0 && return 0  # if the parser get stop token
        return v
    end
end

lgnot!(state::SMARTSParser, qtree::QueryAtom) = lgnot!(state, qtree, atomprop!)
lgnot!(state::SMARTSParser, qtree::QueryBond) = lgnot!(state, qtree, bondsymbol!)
