#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    coords2d,
    coords2d!,
    cartesian


function coords2d(mol::VectorMol, root)
    required_annotation(mol, :Default)
    required_annotation(mol, :Topology)
    zmatrix = [
        -2 nothing nothing nothing nothing nothing nothing;
        -1 -2 1.0 nothing nothing nothing nothing;
        0 -1 1.0 -2 2/3 nothing 1
    ]
    stack = [root]
    done = []
    pred = Dict(root => 0, 0 => -1, -1 => -2)
    ringmap = Dict(0 => Set(), -1 => Set(), -2 => Set())
    cyclemap = Dict(i => c for (i, c) in enumerate(mol.v[:Cycle]))
    merge!(ringmap, cyclemap)
    rings = mol.annotation[:Topology].cycles
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


function cartesian(zmatrix)
    zmatlen = size(zmatrix, 1)
    coords = zeros(Float64, zmatlen, 3)
    coords[1, :] = [0, 0, 0]
    coords[2, :] = [zmatrix[2, 3], 0, 0]
    ang0 = (1 - zmatrix[3, 5]) * pi
    coords[3, :] = vec([sin(ang0) cos(ang0) 0] * zmatrix[3, 3]) + coords[2, :]
    idxmap = Dict(idx => row for (row, idx) in enumerate(zmatrix[:, 1]))
    for i in 4:zmatlen
        p1 = vec3d(coords[idxmap[zmatrix[i, 2]], :])
        p2 = vec3d(coords[idxmap[zmatrix[i, 4]], :])
        p3 = vec3d(coords[idxmap[zmatrix[i, 6]], :])
        v1 = p1 - p2
        v2 = p3 - p2
        v1u = normalize(v1)
        normalv = normalize(cross(v1, v2))
        ang1 = (1 - zmatrix[i, 5]) * pi
        rot1 = rotation(normalv, ang1)
        ang2 = zmatrix[i, 7] * pi
        rot2 = rotation(v1u, ang2)
        len = zmatrix[i, 3]
        coords[i, :] = vec(rot2 * rot1 * v1u) * len + p1
    end
    indexed = hcat(coords, zmatrix[:, 1])
    sorted = sortslices(indexed, dims=1, by=r->r[4])
    filtered = sorted[sorted[:, 4] .> 0, :] # Remove dummy elements
    filtered[:, 1:3]
end
