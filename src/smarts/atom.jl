#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#


# TODO: radical


const SMARTS_ATOM_COND_SYMBOL = Dict(
    'X' => :connectivity,
    'D' => :degree,
    'v' => :valence,
    'H' => :total_hydrogens,
    'r' => :smallest_ring,
    'R' => :ring_count
)

const SMARTS_CHARGE_SIGN = Dict(
    '+' => 1,
    '-' => -1
)


"""
    atomsymbol!(state::SmartsParserState) -> Union{Pair,Nothing}

Atomsymbol <- Br / Cl / [AaBCcNnOoPpSsFI*]
"""
function atomsymbol!(state::Union{SMILESParser,SMARTSParser})
    sym1 = read(state)
    sym2 = lookahead(state, 1)
    if sym1 == 'C' && sym2 == 'l'
        forward!(state, 2)
        return QueryOperator(:and, [
            QueryLiteral(:symbol, :Cl),
            QueryOperator(:not, [QueryLiteral(:isaromatic)])
        ])
    elseif sym1 == 'B' && sym2 == 'r'
        forward!(state, 2)
        return QueryOperator(:and, [
            QueryLiteral(:symbol, :Br),
            QueryOperator(:not, [QueryLiteral(:isaromatic)])
        ])
    elseif sym1 in "BCNOPSFI"
        forward!(state)
        return QueryOperator(:and, [
            QueryLiteral(:symbol, Symbol(sym1)),
            QueryOperator(:not, [QueryLiteral(:isaromatic)])
        ])
    elseif sym1 in "bcnops"
        forward!(state)
        return QueryOperator(:and, [
            QueryLiteral(:symbol, Symbol(uppercase(sym1))),
            QueryLiteral(:isaromatic)
        ])
    elseif state isa SMARTSParser
        if sym1 == 'A'
            forward!(state)
            return QueryOperator(:not, [QueryLiteral(:isaromatic)])
        elseif sym1 == 'a'
            forward!(state)
            return QueryLiteral(:isaromatic)
        elseif sym1 == '*'
            forward!(state)
            return QueryAny(true)
        end
    end
    # return nothing
end



"""
    atomprop!(state::SmartsParserState) -> Union{Pair,Nothing}

AtomProp <- '\$(' RecursiveQuery ')' / Mass / Symbol / AtomNum / Stereo / CHG / [DHRrvX]
"""
function atomprop!(state::Union{SMILESParser,SMARTSParser})
    sym1 = read(state)
    sym2 = lookahead(state, 1)
    if haskey(ATOMSYMBOLMAP, string(uppercase(sym1), sym2))
        # Two letter atoms
        forward!(state, 2)
        if string(uppercase(sym1), sym2) in ("As", "Se")
            isarom = (islowercase(sym1) ? QueryLiteral(:isaromatic)
                : QueryOperator(:not, [QueryLiteral(:isaromatic)]))
            return QueryOperator(:and, [
                QueryLiteral(:symbol, Symbol(uppercase(sym1), sym2)), isarom
            ])
        end
        isuppercase(sym1) || error("aromatic $(sym1) is not supported")
        return QueryLiteral(:symbol, Symbol(sym1, sym2))
    elseif haskey(SMARTS_ATOM_COND_SYMBOL, sym1)
        # Neighbor and ring conditions
        forward!(state)
        if isdigit(sym2)
            num = parse(Int, sym2)
            forward!(state)
        else
            sym1 in ('r', 'R') && return QueryOperator(:not, [QueryLiteral(:ring_count, 0)])
            num = 1
        end
        return QueryLiteral(SMARTS_ATOM_COND_SYMBOL[sym1], num)
    elseif sym1 in ('A', 'a', '*')
        forward!(state)
        sym1 == 'A' && return QueryOperator(:not, [QueryLiteral(:isaromatic)])
        sym1 == 'a' && return QueryLiteral(:isaromatic)
        sym1 == '*' && return QueryAny(true)
    elseif haskey(ATOMSYMBOLMAP, string(uppercase(sym1)))
        # Single letter atoms
        forward!(state)
        if uppercase(sym1) in "BCNOPS"
            isarom = (islowercase(sym1) ? QueryLiteral(:isaromatic)
                : QueryOperator(:not, [QueryLiteral(:isaromatic)]))
            return QueryOperator(:and, [
                QueryLiteral(:symbol, Symbol(uppercase(sym1))), isarom
            ])
        end
        isuppercase(sym1) || error("aromatic $(sym1) is not supported")
        return QueryLiteral(:symbol, Symbol(sym1))
    elseif sym1 == '#'
        # Atomic number
        forward!(state)
        start = state.pos
        while isdigit(lookahead(state, 1))
            forward!(state)
        end
        num = parse(Int, SubString(state.input, start, state.pos))
        forward!(state)
        return QueryLiteral(:symbol, atom_symbol(num))
    elseif sym1 in keys(SMARTS_CHARGE_SIGN)
        # Charge
        forward!(state)
        if isdigit(sym2)
            chg = parse(Int, sym2)
            forward!(state)
        else
            chg = 1
            while read(state) == sym1
                forward!(state)
                chg += 1
            end
        end
        return QueryLiteral(:charge, chg * SMARTS_CHARGE_SIGN[sym1])
    elseif sym1 == '@'
        # Stereo
        # @ -> anticlockwise, @@ -> clockwise, ? -> or not specified
        forward!(state)
        cw = false
        if read(state) == '@'
            forward!(state)
            cw = true
        end
        if read(state) == '?'
            forward!(state)
            return QueryOperator(:not, [QueryLiteral(:stereo, cw ? :anticlockwise : :clockwise)])
        else
            return QueryLiteral(:stereo, cw ? :clockwise : :anticlockwise)
        end
    elseif isdigit(sym1)
        # Isotope
        start = state.pos
        while isdigit(lookahead(state, 1))
            forward!(state)
        end
        num = SubString(state.input, start, state.pos)
        forward!(state)
        return QueryLiteral(:mass, parse(Int, num))
    elseif sym1 == '$' && sym2 == '('
        # Recursive
        forward!(state, 2)
        start = state.pos
        toclose = 1
        while true
            if read(state) == ')'
                if toclose == 1
                    break
                else
                    toclose -= 1
                end
            elseif read(state) == '('
                toclose += 1
            end
            forward!(state)
        end
        q = SubString(state.input, start, state.pos - 1)
        forward!(state)
        return QueryLiteral(:recursive, String(q))
    end
    # return nothing
end


"""
    atom!(state::SmilesParser) -> Vector{SmilesAtom}

Atom <- '[' AtomProp+ ']' / AtomSymbol
"""
function atom!(state::SMILESParser{T,V,E}) where {T,V,E}
    sym1 = read(state)
    if sym1 == '['
        forward!(state)
        q = lghighand!(state, atomprop!)
        symcls = read(state)
        symcls == ']' || error("(atom)! unexpected token: $(symcls)")
        forward!(state)
    else
        q = atomsymbol!(state)
        q === nothing && return V[] # termination token ().\0
    end
    qd = smiles_dict(q)
    default_qd = Dict{Symbol,Any}(
        :symbol => :C,
        :charge => 0,
        :multiplicity => 1,
        :mass => nothing,
        :isaromatic => false,
        :stereo => :unspecified
    )  # default atom
    if haskey(qd, :total_hydrogens)  # explicit hydrogen (e.g. [CH3]) -> hydrogen nodes
        hcnt = qd[:total_hydrogens]
        if !haskey(qd, :symbol)  # special case: [H2]
            qd[:symbol] = :H
            hcnt -= 1
        end
        merge!(default_qd, qd)
        return [V(default_qd), (V(Dict(:symbol => :H)) for _ in 1:hcnt)...]
    else
        merge!(default_qd, qd)
        return [V(default_qd)]
    end
end


"""
    atom!(state::SmartsParser) -> Vector{SmartsAtom}

Atom <- '[' (AtomProp / LogicalOperator)+ ']' / AtomSymbol
"""
function atom!(state::SMARTSParser{T,V,E}) where {T,V,E}
    sym1 = read(state)
    if sym1 == '['
        forward!(state)
        q = lglowand!(state, atomprop!)
        q === nothing && error("(atom!) empty atomprop")
        symcls = read(state)
        symcls == ']' || error("(atom!) unexpected token: $(symcls)")
        forward!(state)
    else
        q = atomsymbol!(state)
        q === nothing && return V[]  # termination token ().\0
    end
    return [V(q)]
end
