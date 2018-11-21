#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    atom!,
    atomquery!,
    atomsymbol!,
    atomprop!


# TODO: radical
# TODO: ring size, ring count


function atom!(state::SmilesParserState)
    """ Atom <- '[' atomprop+ ']' / [0-9] / atomsymbol
    """
    c = read(state)
    if c == '['
        forward!(state)
        a = lghighand!(state, atomprop!)
        ca = read(state)
        @assert ca == ']' "(atom!) unexpected token: $(ca)"
        forward!(state)
        prop = Dict()
        for p in a.second
            if p.first == :and
                for p2 in p.second
                    prop[p2.first] = p2.second
                end
            else
                prop[p.first] = p.second
            end
        end
        return SmilesAtom(
            prop[:Symbol],
            get!(prop, :Charge, 0),
            0,
            get!(prop, :Mass, nothing),
            get!(prop, :Aromatic, false),
            get!(prop, :stereo, nothing)
        )
    elseif isdigit(c)
        num = parse(Int, c)
        forward!(state)
        return :ring => num
    else
        a = atomsymbol!(state)
        if a === nothing
            return
        else
            return SmilesAtom(a[1], 0, 1, nothing, a[2], nothing)
        end
    end
end


function atomquery!(state::SmartsParserState)
    """ Atom <- '[' AtomPropLogicalCond ']' / [0-9] / atomsymbol
    """
    c = read(state)
    if c == '['
        forward!(state)
        q = lglowand!(state, atomprop!)
        cq = read(state)
        @assert cq == ']' "(atomquery!) unexpected token: $(cq)"
        forward!(state)
        return SmartsAtom(q)
    elseif isdigit(c)
        num = parse(Int, c)
        forward!(state)
        return :ring => num
    else
        a = atomsymbol!(state)
        if a === nothing
            return
        elseif a[1] !== nothing && a[2] !== nothing
            q = :and => (:Symbol => a[1], :Aromatic => a[2])
            return SmartsAtom(q)
        elseif a[1] !== nothing
            return SmartsAtom(:Symbol => a[1])
        elseif a[2] !== nothing
            return SmartsAtom(:Aromatic => a[2])
        else
            return SmartsAtom(:Any => true)
        end
    end
end


function atomsymbol!(state::AbstractSmartsParser)
    """ Atomsymbol <- Br / Cl / [AaBCcNnOoPpSsFI*]
    """
    if read(state) == 'C' && lookahead(state, 1) == 'l'
        forward!(state, 2)
        return (:Cl, false)
    elseif read(state) == 'B' && lookahead(state, 1) == 'r'
        forward!(state, 2)
        return (:Br, false)
    elseif read(state) in "BCNOPSFI"
        sym = Symbol(read(state))
        forward!(state)
        return (sym, false)
    elseif read(state) in "cnops"
        sym = Symbol(uppercase(read(state)))
        forward!(state)
        return (sym, true)
    elseif read(state) == 'A'
        forward!(state)
        return (nothing, false)
    elseif read(state) == 'a'
        forward!(state)
        return (nothing, true)
    elseif read(state) == '*'
        forward!(state)
        return (nothing, nothing)
    end
end


const SMARTS_ATOM_COND_SYMBOL = Dict(
    'X' => :Connectivity,
    'D' => :Degree,
    'v' => :Valence,
    'H' => :H_Count
)

const SMARTS_CHARGE_SIGN = Dict(
    '+' => 1,
    '-' => -1
)


function atomprop!(state::AbstractSmartsParser)
    """ AtomReq <- '\$(' RecursiveQuery ')' / ISO / Sym / Num /
        X / v / R / r / Stereo / H / CHG
    """
    c = read(state)
    atomsyms = keys(PERIODIC_TABLE)
    if isuppercase(c)
        # Non-organic atoms
        c2 = lookahead(state, 1)
        c3 = lookahead(state, 2)
        if string(c, c2, c3) in atomsyms
            # Note: three-letter atoms (U-series) are not supported yet
            forward!(state, 3)
            return :Symbol => Symbol(c, c2, c3)
        elseif string(c, c2) in atomsyms
            forward!(state, 2)
            return :Symbol => Symbol(c, c2)
        end
    end
    # Organic atoms
    a = atomsymbol!(state)
    if a !== nothing
        return :and => (:Symbol => a[1], :Aromatic => a[2])
    end
    # Atom properties
    if c in keys(SMARTS_ATOM_COND_SYMBOL)
        forward!(state)
        if isdigit(c2)
            num = parse(Int, c2)
            forward!(state)
        else
            num = 1
        end
        return SMARTS_ATOM_COND_SYMBOL[c] => num
    elseif isuppercase(c) && string(c) in atomsyms
        # Single letter non-organic atoms
        forward!(state, 1)
        return :Symbol => Symbol(c)
    elseif c == '#'
        # Atomic number
        forward!(state)
        start = state.pos
        while isdigit(lookahead(state, 1))
            forward!(state)
        end
        num = parse(Int, SubString(state.input, start, state.pos))
        forward!(state)
        return :Symbol => atomsymbol(num)
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
        return :Charge => chg * SMARTS_CHARGE_SIGN[c]
    elseif c == '@'
        # Stereo
        # @ => 1, @@ => 2, @? => 3, @@? => 4
        s = 1
        if lookahead(state, 1) == '@'
            forward!(state)
            s = 2
        end
        if lookahead(state, 1) == '?'
            forward!(state)
            s += 2
        end
        forward!(state)
        return :stereo => s
    elseif isdigit(c)
        # Isotope
        start = state.pos
        while isdigit(lookahead(state, 1))
            forward!(state)
        end
        num = SubString(state.input, start, state.pos)
        forward!(state)
        return :Mass => parse(Int, num)
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
