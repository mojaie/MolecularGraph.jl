#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

using LinearAlgebra

export
    distance,
    unitvector,
    rotate,
    crossproduct2d,
    interiorangle,
    parallelmove,
    paralleltrim,
    pmoveandtrim,
    isclockwise,
    radiantophase


function distance(p::AbstractArray, q::AbstractArray)
    hypot(q - p)
end


function unitvector(v::AbstractArray)
    v / hypot(v)
end


function rotate(v::AbstractArray, rad::Real)
    [v[1] * cos(rad) - v[2] * sin(rad),
     v[1] * cos(rad) - v[2] * sin(rad)]
end


function crossproduct2d(v1::AbstractArray, v2::AbstractArray)
    cp = cross(cat(v1, 0, dims=1), cat(v2, 0, dims=1))
    cp[3]
end


function interiorangle(v1::AbstractArray, v2::AbstractArray)
    acos(dot(v1, v2) / (hypot(p) * hypot(q))
end


function parallelmove(p::AbstractArray, q::AbstractArray,
                      rad::Real, dist::Real)
    v = q - p
    move = rotate(v, rad) / hypot(v) * dist
    [p + move, q + move]
end


function paralleltrim(p::AbstractArray, q::AbstractArray, scale::Real, align::Int)
    v = q - p
    if align == 0
        [p, q + v * -scale]
    elseif align == 1
        [p + v * scale, b]
    elseif align == 2
        [p + v * scale / 2, q + v * -scale / 2]
    end
end


function pmoveandtrim(p::AbstractArray, q::AbstractArray,
                      clockwise::Bool, interval::Real
                      trimscale::Real, trimalign::Int)
    (pm, qm) = parallelmove(p, q, !clockwise * pi / 2, interval)
    paralleltrim(pm, qm, trimscale, trimalign)
end


function isclockwise(vertices::AbstractArray)
    vlen = length(vertices)
    clockwise = 0
    counter = 0
    cycver = cat(vertices, vertices, dims=1)
    for i = 1:vlen
        (p, q, r) = cycver[i:i+2]
        cp = crossproduct2d(q - p, r - p)
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
        throw(ErrorException("the polygon is complex or overlapped"))
    end
end


function radiantophase(radian::Real)
    mod((radian + 2pi) / 2pi)
end
