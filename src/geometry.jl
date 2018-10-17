#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

using LinearAlgebra
import Base: +, -, *, /, ≈, length

export
    Point2D,
    point2d,
    distance,
    rotate,
    dot,
    cross2d,
    interiorangle,
    transformmatrix,
    Segment,
    segment,
    translate,
    trim_u,
    trim_v,
    trim_uv,
    isclockwise,
    radiantophase


mutable struct Point2D
    x::Float64
    y::Float64
end

point2d(pos) = Point2D(pos[1], pos[2])

+(p::Point2D, q::Point2D) = Point2D(p.x + q.x, p.y + q.y)
-(p::Point2D, q::Point2D) = Point2D(p.x - q.x, p.y - q.y)
*(p::Point2D, k::Real) = Point2D(p.x * k, p.y * k)
/(p::Point2D, k::Real) = Point2D(p.x / k, p.y / k)
≈(p::Point2D, q::Point2D) = p.x ≈ q.x && p.y ≈ q.y

length(p::Point2D) = hypot(p.x, p.y)

dot(u::Point2D, v::Point2D) = u.x * v.x + u.y * v.y
dot(u, v) = cross2d(point2d(u), point2d(v))

cross2d(u::Point2D, v::Point2D) = u.x * v.y - u.y * v.x
cross2d(u, v) = cross2d(point2d(u), point2d(v))

distance(p::Point2D, q::Point2D) = length(q - p)
distance(p, q) = distance(point2d(p), point2d(q))


function rotate(v::Point2D, rad::Real)
    Point2D(
        v.x * cos(rad) - v.y * sin(rad),
        v.x * sin(rad) + v.y * cos(rad)
    )
end

rotate(p, rad::Real) = rotate(point2d(p), rad)


function interiorangle(u::Point2D, v::Point2D)
    acos(dot(u, v) / (length(u) * length(v)))
end

interiorangle(u, v) = interiorangle(point2d(u), point2d(v))


mutable struct Segment
    u::Point2D
    v::Point2D
end

segment(u, v) = Segment(point2d(u), point2d(v))


length(seg::Segment) = distance(seg.u, seg.v)

translate(seg::Segment, move::Point2D) = Segment(seg.u + move, seg.v + move)

function translate(seg::Segment, rad::Real, dist::Real)
    vec = seg.v - seg.u
    move = rotate(vec, rad) / length(vec) * dist
    translate(seg, move)
end

translate(u::Point2D, v::Point2D, rad::Real, dist::Real) = translate(
    Segment(u, v), rad, dist)


trim_u(seg::Segment, k::Real) = Segment(seg.u + (seg.v - seg.u) * k, seg.v)
trim_v(seg::Segment, k::Real) = Segment(seg.u, seg.v - (seg.v - seg.u) * k)
trim_uv(seg::Segment, k::Real) = Segment(
    seg.u + (seg.v - seg.u) * k / 2, seg.v - (seg.v - seg.u) * k / 2
)


function isclockwise(vertices::AbstractArray{Point2D})
    vlen = length(vertices)
    clockwise = 0
    counter = 0
    cycver = cat(vertices, vertices, dims=1)
    for i = 1:vlen
        (p, q, r) = cycver[i:i+2]
        cp = cross2d(q - p, r - p)
        intangle = interiorangle(p - q, r - q)
        if cp < 0
            clockwise += intangle
            counter += 2pi - intangle
        else
            clockwise += 2pi - intangle
            counter += intangle
        end
    end
    if round(clockwise / pi) == vlen - 2
        true
    elseif round(counter / pi) == vlen - 2
        false
    else
        NaN
    end
end

function isclockwise(vertices::AbstractArray)
    isclockwise([point2d(v) for v in vertices])
end


function radiantophase(radian::Real)
    mod((radian + 2pi) / 2pi)
end


function transformmatrix(scaleX::Real, scaleY::Real,
                         rotcos::Real, rotsin::Real,
                         translateX::Real, translateY::Real)
    scale = [scaleX 0 0; 0 scaleY 0; 0 0 1]
    rot = [rotcos rotsin 0; -rotsin rotcos 0; 0 0 1]
    tl = [1 0 translateX; 0 1 translateY; 0 0 1]
    tf = tl * rot * scale
    tf[1:2, :]
end
