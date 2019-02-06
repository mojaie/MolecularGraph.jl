#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    matrix, point, segment,
    x, y, z, u, v, ux, uy, uz, vx, vy, vz,
    coord, x_components, y_components, z_components,
    fmt,
    vector,
    rotationmatrix


import LinearAlgebra: cross
import Formatting: fmt


mutable struct Cartesian3D <: Coordinates
    coords::Matrix{Float64}

    function Cartesian3D(mat)
        size(mat, 2) == 3 || throw(
            DimensionMismatch("Unexpected matrix size $(size(mat))"))
        new(mat)
    end
end


struct Point3D
    coords::Matrix{Float64}
    i::Float64
end


struct Segment3D
    coords::Matrix{Float64}
    u::Int
    v::Int
end


matrix(coords::Cartesian3D) = coords.coords
matrix(point::Point3D) = point.coords[point.i, :]
matrix(segment::Segment3D) = segment.coords[[segment.u, segment.v], :]

point(coords::Cartesian2D, i) = Point3D(matrix(coords), i)
segment(coords::Cartesian2D, u, v) = Segment3D(matrix(coords), u, v)


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


fmt(s::Point3D, expr::String) = (
    fmt(expr, x(s)), fmt(expr, y(s)), fmt(expr, z(s)))

function fmt(s::Segment3D, expr::String)
    return (
        (fmt(expr, ux(s)), fmt(expr, _uy(s)), fmt(expr, _uz(s))),
        (fmt(expr, vx(s)), fmt(expr, _vy(s)), fmt(expr, _vz(s)))
    )
end


_vector(segment::Segment3D) = _v(segment) - _u(segment))
vector(segment::Segment3D) = Point3D(_vector(segment)...)


function rotationmatrix(axis::Cartesian3D, angle)
    (x, y, z) = (x(axis), y(axis), z(axis))
    c = cos(angle)
    s = sin(angle)
    a12 = x * y * (1 - c)
    a13 = x * z * (1 - c)
    a23 = y * z * (1 - c)
    return Matrix{Float64}([
        (c + x^2 * (1 - c)) (a12 - z * s) (a13 + y * s);
        (a12 + z * s) (c + y^2 * (1 - c)) (a23 - x * s);
        (a13 - y * s) (a23 + x * s) (c + z^2 * (1 - c))
    ])
end
