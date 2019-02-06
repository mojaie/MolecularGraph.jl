#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    rawdata, setcoord!, cartesian2d, point, segment, cyclicpath,
    x, y, u, v, ux, uy, vx, vy,
    coord, x_components, y_components,
    vector, interiorangle, translate,
    midpoint, trim_u, trim_v, trim_uv,
    transformmatrix, isclockwise, radiantophase


# import Base: +, -, *


mutable struct Cartesian2D <: Coordinates
    coords::Matrix{Float64}
end


struct Point2D
    coords::Matrix{Float64}
    i::Int
end


struct Segment2D
    coords::Matrix{Float64}
    u::Int
    v::Int
end


struct CyclicPath2D
    coords::Matrix{Float64}
    nodes::Vector{Int}
end


rawdata(coords::Cartesian2D) = coords.coords
rawdata(point::Point2D) = point.coords[point.i, :]
rawdata(segment::Segment2D) = segment.coords[[segment.u, segment.v], :]


function cartesian2d(coords::Matrix{Float64})
    size(coords, 2) == 2 || throw(
        DimensionMismatch("Unexpected matrix size $(size(coords))"))
    Cartesian2D(coords)
end

point(coords::Cartesian2D, i) = Point2D(rawdata(coords), i)
point(segment::Segment2D, i) = Point2D(rawdata(segment), i)
point(x::Float64, y::Float64) = Point2D([x y], 1)
segment(coords::Cartesian2D, u, v) = Segment2D(rawdata(coords), u, v)
segment(coords::Matrix{Float64}) = Segment2D(coords, 1, 2)
cyclicpath(coords::Cartesian2D, nodes) = CyclicPath2D(rawdata(coords), nodes)
cyclicpath(coords::Matrix{Float64}) = CyclicPath2D(coords, 1:size(coords, 1))


x(point::Point2D) = point.coords[point.i, 1]
y(point::Point2D) = point.coords[point.i, 2]

_u(segment::Segment2D) = segment.coords[segment.u, :]
_v(segment::Segment2D) = segment.coords[segment.v, :]
u(segment::Segment2D) = point(segment, 1) # Deprecated
v(segment::Segment2D) = point(segment, 2) # Deprecated
ux(segment::Segment2D) = segment.coords[segment.u, 1]
uy(segment::Segment2D) = segment.coords[segment.u, 2]
vx(segment::Segment2D) = segment.coords[segment.v, 1]
vy(segment::Segment2D) = segment.coords[segment.v, 2]

_coord(coords::Cartesian2D, i::Int) = coords.coords[i, :]
coord(coords::Cartesian2D, i::Int) = point(_coord(coords, i))


function setcoord!(segment::Segment2D, point::Point2D, i::Int)
    segment.coords[i, :] = rawdata(point)
end


x_components(coords::Cartesian2D) = coords.coords[:, 1]
y_components(coords::Cartesian2D) = coords.coords[:, 2]

# +(a::Point2D, b::Point2D) = Point2D(rawdata(a) - rawdata(b), 1)
# -(a::Point2D, b::Point2D) = Point2D(rawdata(a) - rawdata(b), 1)
# *(a::Point2D, f::Float64) = Point2D(rawdata(a) * f, 1)

Formatting.fmt(expr::String, s::Point2D) = (fmt(expr, x(s)), fmt(expr, y(s)))

function Formatting.fmt(expr::String, s::Segment2D)
    return (
        (fmt(expr, ux(s)), fmt(expr, uy(s))),
        (fmt(expr, vx(s)), fmt(expr, vy(s)))
    )
end


_vector(segment::Segment2D) = _v(segment) - _u(segment)
vector(segment::Segment2D) = point(_vector(segment)...)


LinearAlgebra.cross(u::Point2D, v::Point2D) = x(u) * y(v) - y(u) * x(v)
LinearAlgebra.cross(s::Segment2D) = ux(s) * vy(s) - uy(s) * vx(s)

interiorangle(u::Point2D, v::Point2D) = acos(
    dot(rawdata(u), rawdata(v)) / (norm(rawdata(u)) * norm(rawdata(v))))
interiorangle(s::Segment2D) = acos(
    dot(_u(s), _v(s)) / (norm(_u(s)) * norm(_v(s))))


function _translate(s::Segment2D, angle, dist)
    rotation = [cos(angle) -sin(angle); sin(angle) cos(angle)]
    move = rotation * normalize(_v(s) - _u(s)) * dist
    return rawdata(s) .+ transpose(move)
end

translate(s::Segment2D, a, d) = segment(_translate(s, a, d))


_midpoint(s::Segment2D) = (_u(s) + _v(s)) / 2
midpoint(s::Segment2D) = point(_midpoint(s)...)

_trim_u(s::Segment2D, k) = rawdata(s) + [transpose((_v(s) - _u(s)) * k); 0 0]
trim_u(s::Segment2D, k)  = segment(_trim_u(s, k))

_trim_v(s::Segment2D, k) = rawdata(s) + [0 0; transpose((_v(s) - _u(s)) * -k)]
trim_v(s::Segment2D, k)  = segment(_trim_v(s, k))

_trim_uv(s::Segment2D, k) = rawdata(s) + [
    transpose((_v(s) - _u(s)) * k / 2); transpose((_v(s) - _u(s)) * -k / 2)]
trim_uv(s::Segment2D, k)  = segment(_trim_uv(s, k))



function transformmatrix(scale::Point2D, rotate::Point2D, translate::Point2D)
    s = Matrix{Float64}([x(scale) 0 0; 0 y(scale) 0; 0 0 1])
    r = Matrix{Float64}([x(rotate) -y(rotate) 0; y(rotate) x(rotate) 0; 0 0 1])
    t = Matrix{Float64}([1 0 x(translate); 0 1 y(translate); 0 0 1])
    return t * r * s
end


function isclockwise(path::CyclicPath2D)
    vlen = length(path.nodes)
    clockwise = 0
    counter = 0
    cycpath = cat(path.nodes, path.nodes, dims=1)
    for i in 1:vlen
        p = transpose(path.coords[cycpath[i], :])
        q = transpose(path.coords[cycpath[i+1], :])
        r = transpose(path.coords[cycpath[i+2], :])
        cp = cross(segment(vcat(q, r) .- p))
        intangle = interiorangle(segment(vcat(p, r) .- q))
        if cp < 0
            clockwise += intangle
            counter += 2pi - intangle
        else
            clockwise += 2pi - intangle
            counter += intangle
        end
    end
    if round(clockwise / pi) == vlen - 2
        return true
    elseif round(counter / pi) == vlen - 2
        return false
    end
    return # overlapped
end


radiantophase(angle) = mod((angle + 2pi) / 2pi)
