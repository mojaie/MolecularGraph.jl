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
    'r' => :sssrsizes,
    'R' => :sssrcount
)

const SMARTS_CHARGE_SIGN = Dict(
    '+' => 1,
    '-' => -1
)


function atomsymbol!(state::SmartsParserState)
    """ Atomsymbol <- Br / Cl / [AaBCcNnOoPpSsFI*]
    """
    if read(state) == 'C' && lookahead(state, 1) == 'l'
        forward!(state, 2)
        return :and => (:atomsymbol => :Cl, :isaromatic => false)
    elseif read(state) == 'B' && lookahead(state, 1) == 'r'
        forward!(state, 2)
        return :and => (:atomsymbol => :Br, :isaromatic => false)
    elseif read(state) in "BCNOPSFI"
        sym = Symbol(read(state))
        forward!(state)
        return :and => (:atomsymbol => sym, :isaromatic => false)
    elseif read(state) in "cnops"
        sym = Symbol(uppercase(read(state)))
        forward!(state)
        return :and => (:atomsymbol => sym, :isaromatic => true)
    elseif read(state) == 'A'
        forward!(state)
        return :isaromatic => false
    elseif read(state) == 'a'
        forward!(state)
        return :isaromatic => true
    elseif read(state) == '*'
        forward!(state)
        return :any => true
    end
end


function atom!(state::SmilesParser)
    """ Atom <- '[' AtomProp+ ']' / AtomSymbol
    """
    c = read(state)
    if c == '['
        forward!(state)
        a = lghighand!(state, atomprop!)
        ca = read(state)
        @assert ca == ']' "(atom!) unexpected token: $(ca)"
        forward!(state)
        prop = collectand(a)
        @assert haskey(prop, :atomsymbol) || prop[:hydrogenconnected] >= 1
        atoms = [SmilesAtom(
            get(prop, :atomsymbol, :H),
            get(prop, :charge, 0),
            1,
            get(prop, :mass, nothing),
            get(prop, :isaromatic, false),
            get(prop, :stereo, :unspecified)
        )]
        hcnt = get(prop, :hydrogenconnected, 0)
        if !haskey(prop, :atomsymbol)
            hcnt -= 1
        end
        for h in 1:hcnt
            push!(atoms, SmilesAtom(:H))
        end
            
        return atoms
    else
        a = atomsymbol!(state)
        if a === nothing
            return SmilesAtom[]
        else
            sym = a.second[1].second
            arom = a.second[2].second
            return [SmilesAtom(sym, 0, 1, nothing, arom, :unspecified)]
        end
    end
end


function atom!(state::SmartsParser)
    """ Atom <- '[' (AtomProp / LogicalOperator)+ ']' / AtomSymbol
    """
    c = read(state)
    if c == '['
        forward!(state)
        q = lglowand!(state, atomprop!)
        cq = read(state)
        @assert cq == ']' "(atomquery!) unexpected token: $(cq)"
        forward!(state)
        return [SmartsAtom(q)]
    else
        a = atomsymbol!(state)
        if a === nothing
            return SmartsAtom[]
        else
            return [SmartsAtom(a)]
        end
    end
end


function atomprop!(state::SmartsParserState)
    """ AtomProp <- '\$(' RecursiveQuery ')' / Mass / Symbol / AtomNum /
        Stereo / CHG / [DHRrvX]
    """
    c = read(state)
    c2 = lookahead(state, 1)
    if haskey(ATOMSYMBOLMAP, string(uppercase(c), c2))
        # Two letter atoms
        forward!(state, 2)
        if string(uppercase(c), c2) in ("As", "Se")
            return :and => (
                :atomsymbol => Symbol(uppercase(c), c2),
                :isaromatic => islowercase(c)
            )
        end
        @assert isuppercase(c)
        return :atomsymbol => Symbol(c, c2)
    elseif haskey(SMARTS_ATOM_COND_SYMBOL, c)
        # Neighbor and ring conditions
        forward!(state)
        if isdigit(c2)
            num = parse(Int, c2)
            forward!(state)
        else
            c in (:r, :R) && return :not => (:sssrcount => 0)
            num = 1
        end
        return SMARTS_ATOM_COND_SYMBOL[c] => num
    elseif c in ('A', 'a', '*')
        forward!(state)
        c == 'A' && return :isaromatic => false
        c == 'a' && return :isaromatic => true
        c == '*' && return :any => true
    elseif haskey(ATOMSYMBOLMAP, string(uppercase(c)))
        # Single letter atoms
        forward!(state)
        if uppercase(c) in "BCNOPS"
            return :and => (
                :atomsymbol => Symbol(uppercase(c)),
                :isaromatic => islowercase(c)
            )
        end
        @assert isuppercase(c)
        return :atomsymbol => Symbol(c)
    elseif c == '#'
        # Atomic number
        forward!(state)
        start = state.pos
        while isdigit(lookahead(state, 1))
            forward!(state)
        end
        num = parse(Int, SubString(state.input, start, state.pos))
        forward!(state)
        return :atomsymbol => atomsymbol(num)
    elseif c in keys(SMARTS_CHARGE_SIGN)
        # Charge
        forward!(state)
        c2 = read(state)
        if isdigit(c2)
            chg = parse(Int, c2)
            forward!(state)
        else
            chg = 1
            while read(state) == c
                forward!(state)
                chg += 1
            end
        end
        return :charge => chg * SMARTS_CHARGE_SIGN[c]
    elseif c == '@'
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
            return :not => (:stereo => cw ? :anticlockwise : :clockwise)
        else
            return :stereo => cw ? :clockwise : :anticlockwise
        end
    elseif isdigit(c)
        # Isotope
        start = state.pos
        while isdigit(lookahead(state, 1))
            forward!(state)
        end
        num = SubString(state.input, start, state.pos)
        forward!(state)
        return :mass => parse(Int, num)
    elseif c == '$' && lookahead(state, 1) == '('
        # Recursive
        forward!(state, 2)
        start = state.pos
        while read(state) != ')'
            forward!(state)
        end
        q = SubString(state.input, start, state.pos - 1)
        forward!(state)
        return :recursive => q
    end
end
