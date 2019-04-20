#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    Cartesian3D, Segment3D, Point3D,
    rotation


struct Cartesian3D{T<:AbstractMatrix{Float64}} <: Coordinates
    coords::T
end


struct Segment3D{T<:AbstractMatrix{Float64}}
    coords::T
end


struct Point3D{T<:AbstractVector{Float64}}
    pos::T
end



cartesian3d(data::AbstractMatrix{Float64}) = Cartesian3D(data)
cartesian3d(c3d::Cartesian3D, vertices::Vector{Int}
    ) = Cartesian3D(@view c3d.coords[vertices, :])

"""
    cartesian2d(coords::Cartesian3D) -> Cartesian2D

Embed `Cartesian3D` into `Cartesian2D` by just removing z-axis component.
"""
cartesian2d(c3d::Cartesian3D) = cartesian2d(@view c3d.coords[:, 1:2])

segment(c3d::Cartesian3D, u, v) = Segment3D(@view c3d.coords[[u, v], :])
segment(u::Point3D, v::Point3D) = Segment3D(transpose(hcat(u.pos, v.pos)))

point(c3d::Cartesian3D, i::Int) = Point3D(@view c3d.coords[i, :])



x(point::Point3D) = point.pos[1]
y(point::Point3D) = point.pos[2]
z(point::Point3D) = point.pos[3]

u(segment::Segment3D) = Point3D(@view segment.coords[1, :])
v(segment::Segment3D) = Point3D(@view segment.coords[2, :])
ux(segment::Segment3D) = segment.coords[1]
uy(segment::Segment3D) = segment.coords[3]
uz(segment::Segment3D) = segment.coords[5]
vx(segment::Segment3D) = segment.coords[2]
vy(segment::Segment3D) = segment.coords[4]
vz(segment::Segment3D) = segment.coords[6]

x_components(c3d::Cartesian3D) = c3d.coords[:, 1]
y_components(c3d::Cartesian3D) = c3d.coords[:, 2]
z_components(c3d::Cartesian3D) = c3d.coords[:, 3]


vector(segment::Segment3D) = v(segment) - u(segment)


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
