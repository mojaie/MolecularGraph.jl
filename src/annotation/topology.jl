#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    topology!


"""
    topology!(mol::VectorMol)

Assign graph topology parameter vectors to the molecule.
"""
function topology!(mol::VectorMol)
    # Ring membership
    mol[:RingMem] = nodes_cycles(mol)
    # Ring size
    mol[:RingSize] = nodes_cyclesizes(mol)
    # Ring bond or not
    mol[:RingBond] = edges_iscyclemember(mol)
    # Ring bond membership
    mol[:RingBondMem] = edges_cycles(mol)
    # Ring membership count
    mol[:RingMemCount] = nodes_cyclecount(mol)
    # Scaffold membership
    mol[:ScaffoldMem] = two_edge_membership(mol)
    # Component membership
    mol[:ComponentMem] = connected_membership(mol)
end
