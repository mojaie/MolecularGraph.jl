#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    internalcoords


import Formatting: fmt


mutable struct InternalCoords <: Coordinates
    nodes::Matrix{Int} # node1, node2 and node3
    geometry::Matrix{Float64} # distance, angle and dihedral
end


struct PointInternal
    nodes::Matrix{Int}
    geometry::Matrix{Float64}
    i::Int
end


function internalcoords(nodes, geometry)
    size(nodes, 2) == 3 || throw(
        DimensionMismatch("Unexpected matrix size $(size(coords))"))
    size(geometry, 2) == 3 || throw(
        DimensionMismatch("Unexpected matrix size $(size(coords))"))
    InternalCoords(nodes, geometry)
end


rawdata(coords::InternalCoords) = (coords.nodes, coords.geometry)
rawdata(point::PointInternal) = (
    point.nodes[point.i, :], point.geometry[point.i, :])

_point(coords::InternalCoords, i::Int) = cat(nodes[i, :3], geometry[i, :3])
point(coords::InternalCoords, i::Int) = PointInternal(rawdata(coords), i)

node1(coords::InternalCoords, i::Int) = coords.nodes[i, 1]
node2(coords::InternalCoords, i::Int) = coords.nodes[i, 2]
node3(coords::InternalCoords, i::Int) = coords.nodes[i, 3]
distance(coords::InternalCoords, i::Int) = coords.geometry[i, 1]
angle(coords::InternalCoords, i::Int) = coords.geometry[i, 2]
dihedral(coords::InternalCoords, i::Int) = coords.geometry[i, 3]


function to_cartesian(zmatrix)
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
