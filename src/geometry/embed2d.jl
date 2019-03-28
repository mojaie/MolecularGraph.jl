#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    compute2dcoords


"""
    compute2dcoords(mol::MolGraph) -> InternalCoordinates

Compute 2D coordinates of the molecular graph.
"""
function compute2dcoords(mol::MolGraph)
    graph = mol.graph
    fragments = InternalCoordinates[]

    # Extract scaffolds
    scaffoldnodes = biconnected_component(graph) # not 2-edge connected
    scaffolds = nodesubgraph(scaffoldnodes)
    for scaffold in connected_component(scaffolds)
        if is_outerplanar(nodesubgraph(scaffold))
            push!(fragments, outerplaner_embed2d(scaffold))
        else
            push!(fragments, cartesian_embed2d(scaffold))
        end
        # TODO: macrocycle templates
    end

    # Chains
    chainnodes = setdiff(nodeset(graph), scaffoldnodes)
    chains = nodesubgraph(chainnodes)
    for chain in connected_component(chains)
        push!(fragments, chain_embed2d(chain))
    end

    # TODO: Merge fragments


    # TODO: Boundary check and avoid overlap


    return coords
end


"""
    chain_embed2d(graph::UndirectedGraph; kwargs...) -> InternalCoords

Return a 2D embedding of chain(=tree graph).
"""
function chain_embed2d(graph)
    nodes = Set(chains)
    while !isempty(nodes)
        longest = longestpath(chain)
        push!(path, chain(longest))
        setdiff!(nodes, longest)
    end
end


"""
    is_outerplanar(mol::VectorMol) -> Bool

Return whether the graph is outerplanar or not
"""
function is_outerplanar(mol::VectorMol)
    # TODO: generalize and move to graph.planarity
    # The given graph should be a biconnected graph

    # A biconnected braph G is K4 minor if the graphs have 2n-3 and more edges.
    edgecount(mol) > nodecount(mol) * 2 - 3 && return false

    # A biconnected graph G is K2,3 minor if and only if any pair of basic
    # cycles has two or more edges in common.
    seen = Set{Set{Int}}()
    for i findall(mol[bond_ringmem] .> 1)
        m = mol[bond_ringmem][i]
        m in seen && return false
        push!(seen, m)
    end
    return true
end


"""
    outerplaner_embed2d(graph::UndirectedGraph; kwargs...) -> InternalCoords

Return a 2D embedding of the outerplanar graph.

A 2D embedding of an outerplanar graph can be easily determined by Hamiltonian
path traversal.
"""
function outerplaner_embed2d(mol::VecterMol)
    coords = internalcoords()
    root = pop!(nodekeys(mol))
    stack = [root]
    pred = Dict(root => root)
    while !isempty(stack)
        c = pop!(stack)

        p1 = pred[c]
        p2 = pred[p1]
        p3 = pred[p2]
        r = pop!(copy(mol[:RingMem][c]))
        ringsize = pop!(copy(mol[:RingSize][c]))

        for (nbr, bond) in neighbors(mol, c)
            if length(mol[bond_ringmem][bond]) == 1 && !(nbr in keys(pred))
                pred[nbr] = c
            end
        end

        if !isempty(isec3)
            # elongate ring
            angle = 1 - 2 / length(rings[isec3[1]])
            setnodes!(coords, c, [p1, p2, p3], [1,0, angle, 0.0])
        elseif !isempty(isec2)
            # branch ring
            angle = 1 - 2 / length(rings[isec2[1]])
            setnodes!(coords, c, [p1, p2, p3], [1,0, angle, 1.0])
        end
        # TODO: spiro
        nextlevels = []
        pop!(nbrs, pred[c])
        backtracked = length(nextlevels) == 0
        append!(stack, nextlevels)
        push!(done, c)
    end
    coords3d = cartesian(zmatrix)
    coords3d[:, 1:2]
end


"""
    cartesian_embed2d(graph::UndirectedGraph; kwargs...) -> Cartesian2D

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
