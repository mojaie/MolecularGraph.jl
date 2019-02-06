#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export elemental!


function elemental!(mol::VectorMol)
    haskey(mol, :Symbol) && return
    # Symbol
    mol[:Symbol] = getproperty.(nodevector(mol), :symbol)
    # Charge
    mol[:Charge] = getproperty.(nodevector(mol), :charge)
    # Radical
    mol[:Multiplicity] = getproperty.(nodevector(mol), :multiplicity)
    # Bond order
    mol[:BondOrder] = getproperty.(edgevector(mol), :order)

    heavyatoms = zeros(Int, nodecount(mol)) # Number of adjacent heavy atoms
    explhcount = zeros(Int, nodecount(mol)) # Number of explicit hydrogens
    heavyatombonds = zeros(Int, nodecount(mol)) # Total order of adjacent bonds
    for (n, node) in nodesiter(mol)
        for (nbr, e) in neighbors(mol, n)
            if mol[:Symbol][nbr] == :H
                explhcount[n] += 1
            else
                heavyatoms[n] += 1
                heavyatombonds[n] += mol[:BondOrder][e]
            end
        end
    end
    # Degree (graph degree including explicit hydrogens)
    mol[:Degree] = heavyatoms + explhcount
    # Valence
    mol[:Valence] = valence.(mol[:Symbol], mol[:Charge])
    # Number of lone pairs
    mol[:LonePair] = lonepair.(mol[:Symbol], mol[:Charge])
    # Hydrogen count
    hcnt = (v, b) -> v === nothing ? 0 : max(0, v - b)
    mol[:H_Count] = hcnt.(mol[:Valence], heavyatombonds)
    # Connectivity (connection including hydrogens)
    mol[:Connectivity] = heavyatoms + mol[:H_Count]
    # Number of pi electrons
    mol[:Pi] = heavyatombonds - heavyatoms
    # Hydrogen bond donor count
    dc = (sym, h) -> sym in (:N, :O) && h > 0
    mol[:H_Donor] = dc.(mol[:Symbol], mol[:H_Count])
    # Hydrogen bond acceptor count
    ac = (sym, lp) -> lp === nothing ? false : sym in (:N, :O, :F) && lp > 0
    mol[:H_Acceptor] = ac.(mol[:Symbol], mol[:LonePair])
    # Standard molecular weight
    weight = (atom, h) -> atomweight(atom) + H_WEIGHT * h
    implHcount = mol[:H_Count] - explhcount
    mol[:MolWeight] = weight.(nodevector(mol), implHcount)
    return
end


function valence(symbol, charge)
    defs = Dict(
        :H => 1, :B => 3, :C => 4, :N => 3, :O => 2, :F => 1,
        :Si => 4, :P => 3, :S => 2, :Cl => 1,
        :As => 3, :Se => 2, :Br => 1, :I => 1
    )
    num = get(defs, symbol, nothing)
    return num === nothing ? nothing : num + charge
end


function lonepair(symbol, charge)
    defs = Dict(
        :H => 0, :B => -1, :C => 0, :N => 1, :O => 2, :F => 3,
        :Si => 0, :P => 1, :S => 2, :Cl => 3,
        :As => 1, :Se => 2, :Br => 3, :I => 3
    )
    num = get(defs, symbol, nothing)
    return num === nothing ? nothing : num - charge
end
