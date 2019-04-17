#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    rawdata, cartesian3d, point, segment,
    x, y, z, u, v, ux, uy, uz, vx, vy, vz,
    coord, x_components, y_components, z_components,
    vector,
    rotationmatrix


struct Point3D
    coords::Matrix{Float64}
    i::Int
end


struct Segment3D
    coords::Matrix{Float64}
    u::Int
    v::Int
end


rawdata(coords::Cartesian3D) = coords.coords
rawdata(point::Point3D) = point.coords[point.i, :]
rawdata(segment::Segment3D) = segment.coords[[segment.u, segment.v], :]


"""
    cartesian3d(coords::Matrix{Float64}) -> Cartesian3D

Construct a new `Cartesian3D` object from `Array{Float64,2}`. The size of the
array along the second dimension should be 3, or it throws `DimensionMismatch`.
"""
function cartesian3d(coords::Matrix{Float64})
    size(coords, 2) == 3 || throw(
        DimensionMismatch("Unexpected matrix size $(size(coords))"))
    return Cartesian3D(coords)
end

"""
    cartesian3d(coords::InternalCoords) -> Cartesian3D

Embed `InternalCoords` into `Cartesian3D`.
"""
function cartesian3d(coords::InternalCoords)
    rcds = Vector{Float64}[]
    push!(rcds, [0.0, 0.0, 0.0])
    push!(rcds, [distance(coords, 2), 0.0, 0.0])
    ang = pi - angle(coords, 3)
    push!(rcds, [sin(ang), cos(ang), 0.0] * distance(coords, 3) + rcds[2])
    for i in 4:length(coords)
        p1 = rcds[label1(coords, i)]
        p2 = rcds[label2(coords, i)]
        p3 = rcds[label3(coords, i)]
        v1 = p1 - p2
        v2 = p3 - p2
        v1u = normalize(v1)
        normalv = normalize(cross(v1, v2))
        rot1 = rotation(point(normalv...), pi - angle(coords, i))
        rot2 = rotation(point(v1u...), dihedral(coords, i))
        dist = distance(coords, i)
        push!(rcds, rot2 * rot1 * v1u * dist + p1)
    end
    c3d = zeros(Float64, length(coords.nodekeys), 3)
    for (i, rcd) in enumerate(rcds)
        if i in keys(coords.nodekeys)
            c3d[coords.nodekeys[i], :] = rcd
        end
    end
    return Cartesian3D(c3d)
end


_point(coords::Cartesian3D, i::Int) = coords.coords[i, :]
_point(segment::Segment3D, i::Int) = segment.coords[i, :]
point(coords::Cartesian3D, i) = Point3D(rawdata(coords), i)
point(x::Float64, y::Float64, z::Float64) = Point3D([x y z], 1)

segment(coords::Cartesian3D, u, v) = Segment3D(rawdata(coords), u, v)


x(point::Point3D) = point.coords[point.i, 1]
y(point::Point3D) = point.coords[point.i, 2]
z(point::Point3D) = point.coords[point.i, 3]

_u(segment::Segment3D) = segment.coords[segment.u, :]
_v(segment::Segment3D) = segment.coords[segment.v, :]
u(segment::Segment3D) = Point3D(_u(segment)...)
v(segment::Segment3D) = Point3D(_v(segment)...)
ux(segment::Segment3D) = segment.coords[segment.u, 1]
uy(segment::Segment3D) = segment.coords[segment.u, 2]
uz(segment::Segment3D) = segment.coords[segment.u, 3]
vx(segment::Segment3D) = segment.coords[segment.v, 1]
vy(segment::Segment3D) = segment.coords[segment.v, 2]
vz(segment::Segment3D) = segment.coords[segment.v, 3]

_coord(coords::Cartesian3D, i::Int) = coords.coords[i, :]
coord(coords::Cartesian3D, i::Int) = Point3D(_coord(coords, i)...)

x_components(coords::Cartesian3D) = coords.coords[:, 1]
y_components(coords::Cartesian3D) = coords.coords[:, 2]
z_components(coords::Cartesian3D) = coords.coords[:, 3]


_vector(segment::Segment3D) = _v(segment) - _u(segment)
vector(segment::Segment3D) = Point3D(_vector(segment)...)


function rotation(axis::Point3D, angle)
    (x1, y1, z1) = (x(axis), y(axis), z(axis))
    c = cos(angle)
    s = sin(angle)
    a12 = x1 * y1 * (1 - c)
    a13 = x1 * z1 * (1 - c)
    a23 = y1 * z1 * (1 - c)
    return Matrix{Float64}([
        (c + x1^2 * (1 - c)) (a12 - z1 * s) (a13 + y1 * s);
        (a12 + z1 * s) (c + y1^2 * (1 - c)) (a23 - x1 * s);
        (a13 - y1 * s) (a23 + x1 * s) (c + z1^2 * (1 - c))
    ])
end
