#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    vec2d,
    vec3d,

import LinearAlgebra: cross

"""
    coords2d(mol::MolGraph) -> InternalCoordinates

Depth first search based 2D embedding for outerplanar graph.
"""
function coords2d(mol::MolGraph)
    fragments = InternalCoordinates[]

    # Scaffolds
    scaffoldnodes = two_edge_connected(mol)
    scaffolds = nodesubgraph(scaffoldnodes)
    for scaffold in connected_component(scaffolds)
        if isouterplaner(nodesubgraph(scaffold))
            push!(fragments, outerplaner_embed2d(scaffold))
        else
            push!(fragments, cartesian_embed2d(scaffold))
        end
    end

    # Chains
    chainnodes = setdiff(nodekeys(mol), scaffoldnodes)
    chains = nodesubgraph(chainnodes)
    for chain in connected_component(chains)
        push!(fragments, chain_embed2d(chain))
    end

    # Merge fragments


    # Boundary check and avoid overlap

    return coords
end


"""
    chain_embed2d(graph::UDGraph; kwargs...) -> InternalCoords

Return a 2D embedding of chain(=tree graph).
"""
function chain_embed2d(graph)
    nodes = Set(chains)
    while !isempty(nodes)
        pathnodes = longestpath(chain)
    end
end


"""
    outerplaner_embed2d(graph::UDGraph; kwargs...) -> InternalCoords

Return a 2D embedding of the outerplanar graph.

A 2D embedding of an outerplanar graph can be easily determined by depth first
search(DFS) based algorithm.
"""
function outerplaner_embed2d(graph)

end


"""
    cartesian_embed2d(graph::UDGraph; kwargs...) -> Cartesian2D

Cartesian 2D embedding
"""
function cartesian_embed2d(graph)

end


function coords2d(mol::VectorMol, root)
    # Requires :Topology
    zmatrix = [
        -2 nothing nothing nothing nothing nothing nothing;
        -1 -2 1.0 nothing nothing nothing nothing;
        0 -1 1.0 -2 2/3 nothing 1
    ]
    stack = [root]
    done = []
    pred = Dict(root => 0, 0 => -1, -1 => -2)
    ringmap = Dict(0 => Set(), -1 => Set(), -2 => Set())
    cyclemap = Dict(i => c for (i, c) in enumerate(mol[:Cycle]))
    merge!(ringmap, cyclemap)
    rings = mol.annotation[:Topology].rings
    backtracked = false
    while length(stack) > 0
        c = pop!(stack)
        if c in done
            continue
        end
        p1 = pred[c]
        p2 = pred[p1]
        p3 = pred[p2]

        isec3 = collect(intersect(ringmap[c], ringmap[p3]))
        isec2 = collect(intersect(ringmap[c], ringmap[p2]))
        if !isempty(isec3)
            # elongate ring
            angle = 1 - 2 / length(rings[isec3[1]])
            zmatrix = vcat(zmatrix, [c p1 1.0 p2 angle p3 0.0])
        elseif !isempty(isec2)
            # branch ring
            angle = 1 - 2 / length(rings[isec2[1]])
            zmatrix = vcat(zmatrix, [c p1 1.0 p2 angle p3 1.0])
        elseif backtracked
            # branch chain
            zmatrix = vcat(zmatrix, [c p1 1.0 p2 2/3 p3 0.0])
        else
            # elongate chain
            zmatrix = vcat(zmatrix, [c p1 1.0 p2 2/3 p3 1.0])
        end
        # TODO: cis-trans specified
        # TODO: triple bond
        # TODO: spiro
        nextlevels = []
        for n in keys(neighbors(mol, c))
            if n == pred[c]
                continue
            end
            pred[n] = c
            if !(n in done)
                push!(nextlevels, n)
            end
        end
        backtracked = length(nextlevels) == 0
        append!(stack, nextlevels)
        push!(done, c)
    end
    coords3d = cartesian(zmatrix)
    coords3d[:, 1:2]
end

coords2d(mol) = coords2d(mol, 1)
