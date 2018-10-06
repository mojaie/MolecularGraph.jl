#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

using LinearAlgebra
import Base: +, -, *, /, ≈

export
    Point2D,
    distance,
    unitvector,
    rotate,
    dot,
    cross2d,
    interiorangle,
    parallelmove,
    paralleltrim,
    pmoveandtrim,
    isclockwise,
    radiantophase


mutable struct Point2D
    x::Float64
    y::Float64
end

+(p::Point2D, q::Point2D) = Point2D(p.x + q.x, p.y + q.y)
-(p::Point2D, q::Point2D) = Point2D(p.x - q.x, p.y - q.y)
*(p::Point2D, k::Real) = Point2D(p.x * k, p.y * k)
/(p::Point2D, k::Real) = Point2D(p.x / k, p.y / k)
≈(p::Point2D, q::Point2D) = p.x ≈ q.x && p.y ≈ q.y

dot(u::Point2D, v::Point2D) = u.x * v.x + u.y * v.y
dot(u, v) = cross2d(Point2D(u...), Point2D(v...))
cross2d(u::Point2D, v::Point2D) = u.x * v.y - u.y * v.x
cross2d(u, v) = cross2d(Point2D(u...), Point2D(v...))


function distance(p::Point2D, q::Point2D)
    diff = q - p
    hypot(diff.x, diff.y)
end

distance(p, q) = distance(Point2D(p...), Point2D(q...))


function unitvector(v::Point2D)
    v / hypot(v.x, v.y)
end


function rotate(v::Point2D, rad::Real)
    Point2D(
        v.x * cos(rad) - v.y * sin(rad),
        v.x * sin(rad) + v.y * cos(rad)
    )
end

rotate(p, rad::Real) = rotate(Point2D(p...), rad)


function interiorangle(u::Point2D, v::Point2D)
    acos(dot(u, v) / (hypot(u.x, u.y) * hypot(v.x, v.y)))
end

interiorangle(u, v) = interiorangle(Point2D(u...), Point2D(v...))


function parallelmove(p::Point2D, q::Point2D, rad::Real, dist::Real)
    v = q - p
    move = rotate(v, rad) / hypot(v.x, v.y) * dist
    [p + move, q + move]
end

parallelmove(p, q, rad::Real, dist::Real) = parallelmove(
    Point2D(p...), Point2D(q...), rad, dist)


function paralleltrim(p::Point2D, q::Point2D, scale::Real, align::Int=0)
    v = q - p
    if align == 1
        [p, q + v * -scale]
    elseif align == 2
        [p + v * scale, q]
    elseif align == 0
        [p + v * scale / 2, q + v * -scale / 2]
    end
end

paralleltrim(p, q, scale::Real, align::Int=0) = paralleltrim(
    Point2D(p...), Point2D(q...), scale, align)


function pmoveandtrim(p::Point2D, q::Point2D,
                      clockwise::Bool, interval::Real,
                      trimscale::Real, trimalign::Int)
    (pm, qm) = parallelmove(p, q, !clockwise * pi / 2, interval)
    paralleltrim(pm, qm, trimscale, trimalign)
end


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
    isclockwise([Point2D(v...) for v in vertices])
end


function radiantophase(radian::Real)
    mod((radian + 2pi) / 2pi)
end
