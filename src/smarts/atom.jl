#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#


# TODO: radical


const SMARTS_ATOM_COND_SYMBOL = Dict(
    'X' => :connectivity,
    'D' => :nodedegree,
    'v' => :valence,
    'H' => :hydrogenconnected,
    'r' => :smallestsssr,
    'R' => :sssrcount
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
        return (v -> v[1] & ~v[2], [(:symbol, :Cl), (:isaromatic,)])
    elseif sym1 == 'B' && sym2 == 'r'
        forward!(state, 2)
        return (v -> v[1] & ~v[2], [(:symbol, :Br), (:isaromatic,)])
    elseif sym1 in "BCNOPSFI"
        forward!(state)
        return (v -> v[1] & ~v[2], [(:symbol, Symbol(sym1)), (:isaromatic,)])
    elseif sym1 in "bcnops"
        forward!(state)
        return (v -> v[1] & v[2], [(:symbol, Symbol(uppercase(sym1))), (:isaromatic,)])
    elseif sym1 == 'A'
        forward!(state)
        return (v -> ~v[1], [(:isaromatic,)])
    elseif sym1 == 'a'
        forward!(state)
        return (v -> v[1], [(:isaromatic,)])
    elseif sym1 == '*'
        forward!(state)
        return any_query(true)
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
            return (v -> v[1] & v[2], [
                (:symbol, Symbol(uppercase(sym1), sym2)),
                (:isaromatic, islowercase(sym1))
            ])
        end
        isuppercase(sym1) || throw(ErrorException("aromatic $(sym1) is not supported"))
        return (v -> v[1], [(:symbol, Symbol(sym1, sym2))])
    elseif haskey(SMARTS_ATOM_COND_SYMBOL, sym1)
        # Neighbor and ring conditions
        forward!(state)
        if isdigit(sym2)
            num = parse(Int, sym2)
            forward!(state)
        else
            sym1 in ('r', 'R') && return (v -> ~v[1], [(:sssrcount, 0)])
            num = 1
        end
        return (v -> v[1], [(SMARTS_ATOM_COND_SYMBOL[sym1], num)])
    elseif sym1 in ('A', 'a', '*')
        forward!(state)
        sym1 == 'A' && return (v -> ~v[1], [(:isaromatic,)])
        sym1 == 'a' && return (v -> v[1], [(:isaromatic,)])
        sym1 == '*' && return any_query(true)
    elseif haskey(ATOMSYMBOLMAP, string(uppercase(sym1)))
        # Single letter atoms
        forward!(state)
        if uppercase(sym1) in "BCNOPS"
            fml = islowercase(sym1) ? v -> v[1] & v[2] : v -> v[1] & ~v[2]
            return (fml, [(:symbol, Symbol(uppercase(sym1))), (:isaromatic,)])
        end
        isuppercase(sym1) || throw(ErrorException("aromatic $(sym1) is not supported"))
        return (v -> v[1], [(:symbol, Symbol(sym1))])
    elseif sym1 == '#'
        # Atomic number
        forward!(state)
        start = state.pos
        while isdigit(lookahead(state, 1))
            forward!(state)
        end
        num = parse(Int, SubString(state.input, start, state.pos))
        forward!(state)
        return (v -> v[1], [(:symbol, atomsymbol(num))])
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
        return (v -> v[1], [(:charge, chg * SMARTS_CHARGE_SIGN[sym1])])
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
            return (v -> ~v[1], [(:stereo, cw ? :anticlockwise : :clockwise)])
        else
            return (v -> v[1], [(:stereo, cw ? :clockwise : :anticlockwise)])
        end
    elseif isdigit(sym1)
        # Isotope
        start = state.pos
        while isdigit(lookahead(state, 1))
            forward!(state)
        end
        num = SubString(state.input, start, state.pos)
        forward!(state)
        return (v -> v[1], [(:mass, parse(Int, num))])
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
        return (v -> v[1], [(:recursive, String(q))])
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
        symcls == ']' || throw(ErrorException("(atom!) unexpected token: $(symcls)"))
        forward!(state)
    else
        q = atomsymbol!(state)
        q === nothing && return V[] # termination token ().\0
    end
    # :or queries are not included in SMILES, so just accumulate property tuples
    qd = Dict{Symbol,Any}()
    for fml in q[2]  # operator is fixed to :eq so far
        qd[fml[1]] = length(fml) == 1 ? q[1](trues(length(q[2]))) : fml[2]
    end
    if haskey(qd, :hydrogenconnected)  # explicit hydrogen (e.g. [CH3]) -> hydrogen nodes
        hcnt = qd[:hydrogenconnected]
        if !haskey(qd, :symbol)  # special case: [H2]
            qd[:symbol] = :H
            hcnt -= 1
        end
        return [V(qd), (V(Dict(:symbol => :H)) for _ in 1:hcnt)...]
    else
        return [V(qd)]
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
        q === nothing && throw(ErrorException("(atom!) empty atomprop"))
        symcls = read(state)
        symcls == ']' || throw(ErrorException("(atom!) unexpected token: $(symcls)"))
        forward!(state)
    else
        q = atomsymbol!(state)
        q === nothing && return V[]  # termination token ().\0
    end
    return [V(q...)]
end
