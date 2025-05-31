#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

# TODO: Gaussian compatible z-matrix format

struct InternalCoords <: Coordinates
    labels::Matrix{Union{Int,Nothing}} # atom1, atom2 and atom3
    geometry::Matrix{Union{Float64,Nothing}} # distance, angle and dihedral
end

function InternalCoords(size::Int)
    labels = fill(nothing, (size, 3))
    geometry = fill(nothing, (size, 3))
    return InternalCoords(labels, geometry)
end

internalcoords(labels, geometry) = InternalCoords(labels, geometry)


function rotation(axis::Point3d, angle)
    x1, y1, z1 = collect(axis)
    c = cos(angle)
    s = sin(angle)
    a12 = x1 * y1 * (1 - c)
    a13 = x1 * z1 * (1 - c)
    a23 = y1 * z1 * (1 - c)
    return Mat3d([
        (c + x1^2 * (1 - c)) (a12 - z1 * s) (a13 + y1 * s);
        (a12 + z1 * s) (c + y1^2 * (1 - c)) (a23 - x1 * s);
        (a13 - y1 * s) (a23 + x1 * s) (c + z1^2 * (1 - c))
    ])
end


"""
    cartesian3d(coords::InternalCoords) -> Cartesian3D

Embed `InternalCoords` into `Cartesian3D`.
"""
function cartesian3d(cint::InternalCoords, mapping::Dict{Int,Int})
    rcds = Vector{Float64}[]
    push!(rcds, [0.0, 0.0, 0.0])
    push!(rcds, [distance(cint, 2), 0.0, 0.0])
    ang = pi - angle(cint, 3)
    push!(rcds, [sin(ang), cos(ang), 0.0] * distance(cint, 3) + rcds[2])
    for i in 4:length(cint)
        p1 = rcds[label1(cint, i)]
        p2 = rcds[label2(cint, i)]
        p3 = rcds[label3(cint, i)]
        v1 = p1 - p2
        v2 = p3 - p2
        v1u = normalize(v1)
        normalv = normalize(cross(v1, v2))
        rot1 = rotation(point(normalv...), pi - angle(cint, i))
        rot2 = rotation(point(v1u...), dihedral(cint, i))
        dist = distance(cint, i)
        push!(rcds, rot2 * rot1 * v1u * dist + p1)
    end
    c3d = zeros(Float64, length(cint.nodekeys), 3)
    for (i, rcd) in enumerate(rcds)
        if i in keys(cint.nodekeys)
            c3d[cint.nodekeys[i], :] = rcd
        end
    end
    return c3d
end


Base.length(coords::InternalCoords) = size(coords.labels, 1)


label1(coords::InternalCoords, i::Int) = coords.labels[i, 1]
label2(coords::InternalCoords, i::Int) = coords.labels[i, 2]
label3(coords::InternalCoords, i::Int) = coords.labels[i, 3]

distance(coords::InternalCoords, i::Int) = coords.geometry[i, 1]
Base.angle(coords::InternalCoords, i::Int) = coords.geometry[i, 2]
dihedral(coords::InternalCoords, i::Int) = coords.geometry[i, 3]


function setcoord!(coords::InternalCoords, i::Int, labels, geometry)
    coords.labels[i, :] = labels
    coords.geometry[i, :] = geometry
end
