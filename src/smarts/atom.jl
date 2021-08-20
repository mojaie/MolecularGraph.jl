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
function atomsymbol!(state::SmartsParserState)
    sym1 = read(state)
    sym2 = lookahead(state, 1)
    if sym1 == 'C' && sym2 == 'l'
        forward!(state, 2)
        return QueryFormula(
            :and, Set([QueryFormula(:atomsymbol, :Cl), QueryFormula(:isaromatic, false)]))
    elseif sym1 == 'B' && sym2 == 'r'
        forward!(state, 2)
        return QueryFormula(
            :and, Set([QueryFormula(:atomsymbol, :Br), QueryFormula(:isaromatic, false)]))
    elseif sym1 in "BCNOPSFI"
        forward!(state)
        return QueryFormula(
            :and, Set([QueryFormula(:atomsymbol, Symbol(sym1)), QueryFormula(:isaromatic, false)]))
    elseif sym1 in "bcnops"
        forward!(state)
        return QueryFormula(
            :and, Set([
                QueryFormula(:atomsymbol, Symbol(uppercase(sym1))),
                QueryFormula(:isaromatic, true)
            ]))
    elseif sym1 == 'A'
        forward!(state)
        return QueryFormula(:isaromatic, false)
    elseif sym1 == 'a'
        forward!(state)
        return QueryFormula(:isaromatic, true)
    elseif sym1 == '*'
        forward!(state)
        return QueryFormula(:any, true)
    end
    # return nothing
end



"""
    atomprop!(state::SmartsParserState) -> Union{Pair,Nothing}

AtomProp <- '\$(' RecursiveQuery ')' / Mass / Symbol / AtomNum / Stereo / CHG / [DHRrvX]
"""
function atomprop!(state::SmartsParserState)
    sym1 = read(state)
    sym2 = lookahead(state, 1)
    if haskey(ATOMSYMBOLMAP, string(uppercase(sym1), sym2))
        # Two letter atoms
        forward!(state, 2)
        if string(uppercase(sym1), sym2) in ("As", "Se")
            return QueryFormula(
                :and, Set([
                    QueryFormula(:atomsymbol, Symbol(uppercase(sym1), sym2)),
                    QueryFormula(:isaromatic, islowercase(sym1))
                ]))
        end
        @assert isuppercase(sym1)
        return QueryFormula(:atomsymbol, Symbol(sym1, sym2))
    elseif haskey(SMARTS_ATOM_COND_SYMBOL, sym1)
        # Neighbor and ring conditions
        forward!(state)
        if isdigit(sym2)
            num = parse(Int, sym2)
            forward!(state)
        else
            sym1 in ('r', 'R') && return QueryFormula(:not, QueryFormula(:sssrcount, 0))
            num = 1
        end
        return QueryFormula(SMARTS_ATOM_COND_SYMBOL[sym1], num)
    elseif sym1 in ('A', 'a', '*')
        forward!(state)
        sym1 == 'A' && return QueryFormula(:isaromatic, false)
        sym1 == 'a' && return QueryFormula(:isaromatic, true)
        sym1 == '*' && return QueryFormula(:any, true)
    elseif haskey(ATOMSYMBOLMAP, string(uppercase(sym1)))
        # Single letter atoms
        forward!(state)
        if uppercase(sym1) in "BCNOPS"
            return QueryFormula(
                :and, Set([
                    QueryFormula(:atomsymbol, Symbol(uppercase(sym1))),
                    QueryFormula(:isaromatic, islowercase(sym1))
                ]))
        end
        @assert isuppercase(sym1)
        return QueryFormula(:atomsymbol, Symbol(sym1))
    elseif sym1 == '#'
        # Atomic number
        forward!(state)
        start = state.pos
        while isdigit(lookahead(state, 1))
            forward!(state)
        end
        num = parse(Int, SubString(state.input, start, state.pos))
        forward!(state)
        return QueryFormula(:atomsymbol, atomsymbol(num))
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
        return QueryFormula(:charge, chg * SMARTS_CHARGE_SIGN[sym1])
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
            return QueryFormula(:not, QueryFormula(:stereo, cw ? :anticlockwise : :clockwise))
        else
            return QueryFormula(:stereo, cw ? :clockwise : :anticlockwise)
        end
    elseif isdigit(sym1)
        # Isotope
        start = state.pos
        while isdigit(lookahead(state, 1))
            forward!(state)
        end
        num = SubString(state.input, start, state.pos)
        forward!(state)
        return QueryFormula(:mass, parse(Int, num))
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
        return QueryFormula(:recursive, q)
    end
    # return nothing
end



"""
    atom!(state::SmilesParser) -> Vector{SmilesAtom}

Atom <- '[' AtomProp+ ']' / AtomSymbol
"""
function atom!(state::SmilesParser)
    sym1 = read(state)
    if sym1 == '['
        forward!(state)
        fml = lghighand!(state, atomprop!)
        fml = tidyformula(fml)
        symcls = read(state)
        @assert symcls == ']' "(atom!) unexpected token: $(symcls)"
        forward!(state)
        @assert findformula(fml, :atomsymbol) !== nothing || findformula(fml, :hydrogenconnected) >= 1
        atoms = [SmilesAtom(
            something(findformula(fml, :atomsymbol), :H),
            something(findformula(fml, :charge), 0),
            1,  # TODO: multiplicity
            findformula(fml, :mass),
            something(findformula(fml, :isaromatic), false),
            something(findformula(fml, :stereo), :unspecified)
        )]
        hcnt = something(findformula(fml, :hydrogenconnected), 0)
        if findformula(fml, :atomsymbol) === nothing
            hcnt -= 1
        end
        for h in 1:hcnt
            push!(atoms, SmilesAtom(:H))
        end
        return atoms
    else
        fml = atomsymbol!(state)
        if fml === nothing
            return SmilesAtom[]
        else
            asym = findformula(fml, :atomsymbol)
            arom = findformula(fml, :isaromatic)
            return [SmilesAtom(asym, 0, 1, nothing, arom, :unspecified)]
        end
    end
end


"""
    atom!(state::SmartsParser) -> Vector{SmartsAtom}

Atom <- '[' (AtomProp / LogicalOperator)+ ']' / AtomSymbol
"""
function atom!(state::SmartsParser)
    sym1 = read(state)
    if sym1 == '['
        forward!(state)
        fml = lglowand!(state, atomprop!)
        @assert fml !== nothing "(atom!) empty atomprop"
        fml = tidyformula(fml)
        symcls = read(state)
        @assert symcls == ']' "(atom!) unexpected token: $(symcls)"
        forward!(state)
        return [SmartsAtom(fml)]
    else
        fml = atomsymbol!(state)
        if fml === nothing
            return SmartsAtom[]
        else
            return [SmartsAtom(fml)]
        end
    end
end
