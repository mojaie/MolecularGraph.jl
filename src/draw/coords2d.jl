#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    coords2d,
    coords2d!,
    cartesian

using LinearAlgebra


function coords2d(mol::MolecularGraph, root::Integer)
    required_descriptor(mol, "Valence")
    required_descriptor(mol, "Topology")
    zmatrix = [
        -2 nothing nothing nothing nothing nothing nothing;
        -1 -2 1.0 nothing nothing nothing nothing;
        0 -1 1.0 -2 2/3 nothing 1
    ]
    stack = [root]
    done = []
    pred = Dict(root => 0, 0 => -1, -1 => -2)
    ringmap = Dict{Int64, Array}(a.index => Int64[] for a in atomvector(mol))
    merge!(ringmap, Dict(0 => Int64[], -1 => Int64[], -2 => Int64[]))
    for (i, ring) in enumerate(mol.rings)
        for r in ring.arr
            push!(ringmap[r], i)
        end
    end
    backtracked = false
    while length(stack) > 0
        c = pop!(stack)
        if c in done
            continue
        end
        p1 = pred[c]
        p2 = pred[p1]
        p3 = pred[p2]

        isec3 = intersect(ringmap[c], ringmap[p3])
        isec2 = intersect(ringmap[c], ringmap[p2])
        if !isempty(isec3)
            # elongate ring
            ringsize = length(mol.rings[isec3[1]])
            angle = (ringsize - 2) / ringsize
            zmatrix = vcat(zmatrix, [c p1 1.0 p2 angle p3 0.0])
        elseif !isempty(isec2)
            # branch ring
            ringsize = length(mol.rings[isec2[1]])
            angle = (ringsize - 2) / ringsize
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
    zmatrix
end

coords2d(mol) = coords2d(mol, 1)


function coords2d!(mol::MolecularGraph)
    zmatrix = coords2d(mol)
    coords3d = cartesian(zmatrix)
    for (i, atom) in enumerate(atomvector(mol))
        atom.coords = tuple(coords3d[i, 1:2]...)
    end
end


function cartesian(zmatrix::AbstractArray)
    zmatlen = size(zmatrix, 1)
    coords = zeros(Float32, zmatlen, 3)
    coords[1, :] = [0, 0, 0]
    coords[2, :] = [zmatrix[2, 3], 0, 0]
    ang0 = (1 - zmatrix[3, 5]) * pi
    coords[3, :] = vec([sin(ang0) cos(ang0) 0] .* zmatrix[3, 3]) .+ coords[2, :]
    idxmap = Dict(idx => row for (row, idx) in enumerate(zmatrix[:, 1]))
    for i in 4:zmatlen
        p1 = coords[idxmap[zmatrix[i, 2]], :]
        p2 = coords[idxmap[zmatrix[i, 4]], :]
        p3 = coords[idxmap[zmatrix[i, 6]], :]
        v1 = p1 .- p2
        v2 = p3 .- p2
        v1u = v1 / hypot(v1...)
        normal = cross(v1, v2)
        normalu = normal / hypot(normal...)
        ang1 = (1 - zmatrix[i, 5]) * pi
        rot1 = rotationmatrix(normalu, ang1)
        ang2 = zmatrix[i, 7] * pi
        rot2 = rotationmatrix(v1u, ang2)
        len = zmatrix[i, 3]
        coords[i, :] = vec(rot2 * rot1 * reshape(v1u, (3, 1))) .* len .+ p1
    end
    indexed = hcat(coords, zmatrix[:, 1])
    sorted = sortslices(indexed, dims=1, by=r->r[4])
    filtered = sorted[sorted[:, 4] .> 0, :] # Remove dummy elements
    filtered[:, 1:3]
end
