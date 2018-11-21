#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    smilestomol,
    parsesmiles,
    tokenize,
    parsesmilesgroup!,
    parsesmilestoken,
    parsesmilesatom,
    parsesmilesbond


function smilestomol(smiles)
    mol = parsesmiles(smiles)
    vmol = vectormol(mol)
    default_annotation!(vmol)
    vmol
end


function parsesmiles(smiles)
    mol = SMILES()
    tokens = tokenize(smiles)
    ringclose = Dict()
    parsesmilesgroup!(mol, ringclose, tokens, 0) # Recursive
    return mol
end


function tokenize(smiles)
    tokens = []
    i = 1
    while i != lastindex(smiles) + 1
        m = match(
            r"^[\.\(\)]|^([=#:/\\]?(\[.+?\]|Br|Cl|[BCcNnOoPpSsFI])([=#:/\\]?[1-9]+)*)",
            smiles[i:end]
        )
        if m === nothing
            throw(MolParseError("Invalid SMILES"))
        end
        i += length(m.match)
        push!(tokens, m.match)
    end
    tokens
end


function parsesmilesgroup!(mol, ringclose, tokens, base)
    pos = base
    nobond = base == 0
    while length(tokens) > 0
        token = popfirst!(tokens)
        if token == "("
            parsesmilesgroup!(mol, ringclose, tokens, pos)
            continue
        elseif token == ")"
            break
        elseif token == "."
            nobond = true
            continue
        end

        newpos = length(mol.graph.nodes) + 1
        # Atom
        # println(objectid(mol), " ", pos, " ", newpos, " ", token)
        (atok, btok, ringc) = parsesmilestoken(token)
        atom = SmilesAtom(parsesmilesatom(atok)...)
        updateatom!(mol, atom, newpos)

        # Ring bond
        if ringc !== nothing
            rtoks = []
            i = 1
            while i != lastindex(ringc) + 1
                m = match(r"^[=#:/\\]?[1-9]", ringc[i:end])
                i += length(m.match)
                push!(rtoks, m.match)
            end
            for rtok in rtoks
                num = rtok[end]
                bondtype = length(rtok) == 2 ? SubString(rtok, 1, 1) : nothing
                if num in keys(ringclose)
                    (upos, ubd) = ringclose[num]
                    rtype = ubd === nothing ? bondtype : ubd
                    bond = SmilesBond(upos, newpos, parsesmilesbond(rtype)...)
                    updatebond!(mol, bond, length(mol.graph.edges) + 1)
                    delete!(ringclose, num)
                else
                    ringclose[num] = (newpos, bondtype)
                end
            end
        end

        # Bond
        if nobond
            nobond = false
        else
            bond = SmilesBond(pos, newpos, parsesmilesbond(btok)...)
            updatebond!(mol, bond, length(mol.graph.edges) + 1)
        end
        pos = newpos
    end
    return
end


function parsesmilestoken(token)
    m = match(
        r"^([=#:/\\])?(\[.+?\]|Br|Cl|[BCcNnOoPpSsFI])([1-9=#:/\\]*)",
        token
    )
    bond = m.captures[1]
    atom = replace(m.captures[2], r"\[(.+?)\]" => s"\1")
    ringc = m.captures[3] == "" ? nothing : m.captures[3]
    (atom, bond, ringc)
end


const SMILES_CHARGE_TABLE = Dict(
    nothing => 0, "+" => 1, "++" => 2, "+++" => 3, "++++" => 4,
    "-" => -1, "--" => -2, "---" => -3, "----" => -4,
    "+" => 1, "+2" => 2, "+3" => 3, "+4" => 4,
    "-" => -1, "-2" => -2, "-3" => -3, "-4" => -4
)
const SMILES_H_COUNT_TABLE = Dict(
    nothing => 0, "H" => 1, "H2" => 2, "H3" => 3, "H4" => 4)

const SMILES_STEREO_TABLE = Dict("" => nothing, "@" => 1, "@@" => 2)


function parsesmilesatom(token)
    m = match(
        r"^([0-9]*)([cnposA-Z][a-z]*)(@*)(H[2-4]?)?([\+\-]+[1-4]?)?",
        token
    )
    sym = Symbol(uppercasefirst(m.captures[2]))
    charge = SMILES_CHARGE_TABLE[m.captures[5]]
    multi = 1 # TODO: OpenBabel radical notation
    mass = m.captures[1] == "" ? nothing : parse(Float64, m.captures[1])
    stereo = SMILES_STEREO_TABLE[m.captures[3]]
    aromatic = islowercase(m.captures[2][1])
    # h_count will be annotated later
    # hcount = SMILES_H_COUNT_TABLE[m.captures[4]])
    [sym, charge, multi, mass, aromatic, stereo]
end


const SMILES_ORDER_TABLE = Dict(
    nothing => 1, "=" => 2, "#" => 3, ":" => 1, "/" => 1, "\\" => 1
)
const SMILES_CISTRANS_TABLE = Dict("/" => 1, "\\" => 2)


function parsesmilesbond(token)
    order = SMILES_ORDER_TABLE[token]
    aromatic = token == ":"
    cistrans = get(SMILES_CISTRANS_TABLE, token, 0)
    [order, aromatic, cistrans]
end
