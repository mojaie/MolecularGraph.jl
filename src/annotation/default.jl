#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export default_annot!


struct DefaultAnnot <: Annotation end


function default_annot!(mol::VectorMol)
    # Symbol
    mol.v[:Symbol] = [atom.symbol for atom in mol.graph.nodes]
    # Charge
    mol.v[:Charge] = [atom.charge for atom in mol.graph.nodes]
    # Radical
    mol.v[:Multiplicity] = [atom.multiplicity for atom in mol.graph.nodes]
    # Bond order
    mol.v[:BondOrder] = [bond.order for bond in mol.graph.edges]

    nbrords = Vector{Int}[]
    for (i, nbr) in enumerate(mol.graph.adjacency)
        push!(nbrords, [mol.v[:BondOrder][b] for b in values(nbr)])
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
    mol.v[:MolWeight] = molweight.(mol.graph.nodes, mol.v[:H_Count])
    mol.annotation[:Default] = DefaultAnnot()
    return
end


function lonepair(symbol, charge)
    defs = Dict(
        :B => -1, :C => 0, :N => 1, :O => 2, :F => 3,
        :Si => 0, :P => 1, :S => 2, :Cl => 3,
        :As => 1, :Se => 2, :Br => 3, :I => 3
    )
    num = get(defs, symbol, nothing)
    num === nothing ? nothing : num - charge
end


function h_count(valence, lonepair)
    lonepair === nothing ? 0 : max(0, 4 - abs(lonepair) - valence)
end


function h_donor(symbol, h_count)
    symbol in (:N, :O) && h_count > 0
end


function h_acceptor(symbol, lonepair)
    symbol in (:N, :O, :F) && lonepair > 0
end


function molweight(atom, h_count)
    atomweight(atom) + H_WEIGHT * h_count
end
