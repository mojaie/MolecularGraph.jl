#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export elemental!


struct Elemental <: Annotation end


# TODO: consider SMILES aromatic

function elemental!(mol::VectorMol)
    # Symbol
    mol.v[:Symbol] = [atom.symbol for (i, atom) in nodesiter(mol)]
    # Charge
    mol.v[:Charge] = [atom.charge for (i, atom) in nodesiter(mol)]
    # Radical
    mol.v[:Multiplicity] = [atom.multiplicity for (i, atom) in nodesiter(mol)]
    # Bond order
    mol.v[:BondOrder] = [bond.order for (i, bond) in edgesiter(mol)]

    nbrords = Vector{Int}[]
    for (n, node) in nodesiter(mol)
        orders = Int[]
        for (nbr, e) in neighbors(mol, n)
            if mol.v[:Symbol][nbr] != :H
                push!(orders, mol.v[:BondOrder][e])
            end
        end
        push!(nbrords, orders)
    end
    # Degree (connection without hydrogen)
    mol.v[:Degree] = length.(nbrords)
    # Valence
    sm = arr -> reduce(+, arr; init=0)  # TODO: sum of empty arr is invalid
    mol.v[:Valence] = sm.(nbrords)
    # Number of pi electrons
    mol.v[:Pi] = mol.v[:Valence] - mol.v[:Degree]
    # Number of lone pairs
    mol.v[:LonePair] = lonepair.(mol.v[:Symbol], mol.v[:Charge])
    # Hydrogen count
    mol.v[:H_Count] = h_count.(mol.v[:Valence], mol.v[:LonePair])
    # Connectivity (connection including hydrogen)
    mol.v[:Connectivity] = mol.v[:Degree] + mol.v[:H_Count]
    # Hydrogen bond donor count
    mol.v[:H_Donor] = h_donor.(mol.v[:Symbol], mol.v[:H_Count])
    # Hydrogen bond acceptor count
    mol.v[:H_Acceptor] = h_acceptor.(mol.v[:Symbol], mol.v[:LonePair])
    # Molecular weight including neighbor hydrogens
    mol.v[:MolWeight] = atomhweight.(
        [atom for (i, atom) in nodesiter(mol)], mol.v[:H_Count])
    mol.annotation[:Elemental] = Elemental()
    return
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


function h_count(valence, lonepair)
    return lonepair === nothing ? 0 : max(0, 4 - abs(lonepair) - valence)
end


function h_donor(symbol, h_count)
    return symbol in (:N, :O) && h_count > 0
end


function h_acceptor(symbol, lonepair)
    return symbol in (:N, :O, :F) && lonepair > 0
end


function atomhweight(atom, h_count)
    if atom.symbol == :H
        return 0  # avoid double count
    else
        return atomweight(atom) + H_WEIGHT * h_count
    end
end
