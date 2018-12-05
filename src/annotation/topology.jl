#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    MolGraphTopology,
    molgraph_topology!


struct MolGraphTopology <: Annotation
    rings::Vector{Vector{Int}}
    scaffolds::Vector{Set{Int}}
    molecules::Vector{Set{Int}}
end


function molgraph_topology!(mol::VectorMol)
    # TODO: shortcut function
    # minimumcycles is derived from 2-edge connected components
    # 2-edge connected components is derived from connected components
    components = connected_components(mol.graph)  # Vector{Set{Int}}
    scaffolds = two_edge_connected(mol.graph)  # Vector{Set{Int}}
    rings = minimumcycles(mol.graph)  # Vector{Vector{Int}}

    mol.annotation[:Topology] = MolGraphTopology(rings, scaffolds, components)
    # Ring membership
    mol.v[:RingMem] = Set{Int}[Set() for i in 1:atomcount(mol)] # dont use fill
    # Ring size
    mol.v[:RingSize] = Set{Int}[Set() for i in 1:atomcount(mol)]
    # Ring bond or not
    mol.v[:RingBond] = falses(bondcount(mol))
    for (i, ring) in enumerate(rings)
        size = length(ring)
        for n in ring
            push!(mol.v[:RingMem][n], i)
            push!(mol.v[:RingSize][n], size)
        end
        sub = inducedsubgraph(mol.graph, Set(ring))
        for e in edgekeys(sub)
            mol.v[:RingBond][e] = true
        end
    end
    # Ring membership count
    mol.v[:RingMemCount] = length.(mol.v[:RingMem])

    # Scaffold membership
    mol.v[:ScaffoldMem] = Vector{Union{Nothing,Int}}(nothing, atomcount(mol))
    for (i, scaffold) in enumerate(scaffolds)
        for s in scaffold
            mol.v[:ScaffoldMem][s] = i
        end
    end
    # Component membership
    mol.v[:ComponentMem] = Vector{Union{Nothing,Int}}(nothing, atomcount(mol))
    for (i, component) in enumerate(components)
        for c in component
            mol.v[:ComponentMem][c] = i
        end
    end
end
