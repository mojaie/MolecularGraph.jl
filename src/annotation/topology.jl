#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    Topology,
    topology!


struct Topology <: Annotation
    rings::Vector{Vector{Int}}
    scaffolds::Vector{Vector{Int}}
    molecules::Vector{Vector{Int}}
end


"""
    topology!(mol::VectorMol)

Assign graph topology parameter vectors to the molecule.
"""
function topology!(mol::VectorMol)
    haskey(mol.annotation, :Topology) && return
    # TODO: shortcut function
    # minimumcycles is derived from 2-edge connected components
    # 2-edge connected components is derived from connected components
    components = connected_components(mol)  # Vector{Set{Int}}
    scaffolds = two_edge_connected(mol)  # Vector{Set{Int}}
    rings = mincycles(mol)  # Vector{Vector{Int}}

    mol.annotation[:Topology] = Topology(rings, scaffolds, components)
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
