#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#


# TODO: migrate to GeometryBasics

export
    Point2D, Point3D, Segment,
    toarray, x_components, y_components, z_components,
    distance, unitvector, midpoint,
    translate, trim_u, trim_v, trim_uv,
    cross2d, interiorangle, isclockwise,
    transformmatrix, rotation


struct Point2D <: MGPoint
    x::Float64
    y::Float64
end

struct Point3D <: MGPoint
    x::Float64
    y::Float64
    z::Float64
end

struct Segment{T<:MGPoint}
    u::T
    v::T
end


Point2D(x::Real, y::Real) = Point2D(float(x), float(y))
Point2D(t::Tuple{T,T}) where {T<:Real} = Point2D(float(t[1]), float(t[2]))
Point2D(arr::AbstractArray{T}
    ) where {T<:Real} = Point2D(float(arr[1]), float(arr[2]))
Point2D(coords::AbstractMatrix{T}, i::Int
    ) where {T<:Real} = Point2D(coords[i, 1:2])

Point3D(x::Real, y::Real, z::Real) = Point3D(float(x), float(y), float(z))
Point3D(t::Tuple{T,T,T}
    ) where {T<:Real} = Point3D(float(t[1]), float(t[2]), float(t[3]))
Point3D(arr::AbstractArray{T}
    ) where {T<:Real} = Point3D(float(arr[1]), float(arr[2]), float(arr[3]))
Point3D(coords::AbstractMatrix{T}, i::Int
    ) where {T<:Real} = Point3D(coords[i, 1:3])

Segment{T}(coords::AbstractMatrix{Float64}, u::Int, v::Int
    ) where {T<:MGPoint} = Segment(T(coords, u), T(coords, v))


Base.:+(u::Point2D, v::Point2D) = Point2D(u.x + v.x,  u.y + v.y)
Base.:+(p::Point2D, transl) = Point2D(p.x + transl[1],  p.y + transl[2])
Base.:+(u::Point3D, v::Point3D) = Point3D(u.x + v.x,  u.y + v.y,  u.z + v.z)
Base.:+(p::Point3D, transl
    ) = Point3D(p.x + transl[1],  p.y + transl[2],  p.z + transl[3])
Base.:+(s::Segment, transl) = Segment(s.u + transl, s.v + transl)

Base.:-(u::Point2D, v::Point2D) = Point2D(u.x - v.x,  u.y - v.y)
Base.:-(p::Point2D, transl) = Point2D(p.x - transl[1],  p.y - transl[2])
Base.:-(u::Point3D, v::Point3D) = Point3D(u.x - v.x,  u.y - v.y,  u.z - v.z)
Base.:-(p::Point3D, transl
    ) = Point3D(p.x - transl[1],  p.y - transl[2],  p.y - transl[3])
Base.:-(s::Segment, transl) = Segment(s.u - transl, s.v - transl)

Base.:*(p::Point2D, f::Float64) = Point2D(p.x * f, p.y * f)
Base.:*(f::Float64, p::Point2D) = Point2D(p.x * f, p.y * f)
Base.:*(p::Point3D, f::Float64) = Point3D(p.x * f, p.y * f, p.z * f)
Base.:*(f::Float64, p::Point3D) = Point3D(p.x * f, p.y * f, p.z * f)
Base.:*(s::Segment, f::Real) = Segment(s.u * f, s.v * f)
Base.:*(f::Real, s::Segment) = Segment(s.u * f, s.v * f)

toarray(p::Point2D) = [p.x, p.y]
toarray(p::Point3D) = [p.x, p.y, p.z]
toarray(u::MGPoint, v::MGPoint) = transpose(hcat(toarray(u), toarray(v)))
toarray(s::Segment) = toarray(s.u, s.v)
toarray(coords::AbstractMatrix{T}) where {T<:Real} = coords
toarray(coords::AbstractMatrix{T}, i::Int
    ) where {T<:Real} = @view coords[i, :]
toarray(coords::AbstractMatrix{T}, u::Int, v::Int
    ) where {T<:Real} = @view coords[[u, v], :]
toarray(coords::AbstractMatrix{T}, vertices::Vector{Int}
    ) where {T<:Real} = @view coords[vertices, :]

x_components(coords::AbstractMatrix{T}) where {T<:Real} = @view coords[:, 1]
y_components(coords::AbstractMatrix{T}) where {T<:Real} = @view coords[:, 2]
z_components(coords::AbstractMatrix{T}) where {T<:Real} = @view coords[:, 3]


LinearAlgebra.norm(p::MGPoint) = norm(toarray(p))
LinearAlgebra.normalize(p::MGPoint) = Point2D(normalize(toarray(p)))
LinearAlgebra.dot(u::MGPoint, v::MGPoint) = dot(toarray(u), toarray(v))


"""
    distance(u::MGPoint, v::MGPoint) -> Float64
    distance(s::Segment) -> Float64
    distance(coords::AbstractMatrix{T}) where {T<:Real} -> Float64

Return distance between two endpoints.
"""
distance(u::MGPoint, v::MGPoint) = norm(v - u)
distance(s::Segment) = distance(s.u, s.v)
distance(coords::AbstractMatrix{T}
    ) where {T<:Real} = norm(coords[2, :] - coords[1, :])


"""
    unitvector(u::MGPoint, v::MGPoint) -> MGPoint
    unitvector(s::Segment) -> MGPoint
    unitvector(coords::AbstractMatrix{T}) where {T<:Real} -> AbstractMatrix

Return u -> v vector of length 1.
"""
unitvector(u::MGPoint, v::MGPoint) = normalize(v - u)
unitvector(s::Segment) = unitvector(s.u, s.v)
unitvector(coords::AbstractMatrix{T}
    ) where {T<:Real} = normalize(coords[2, :] - coords[1, :])


"""
    midpoint(u::MGPoint, v::MGPoint) -> MGPoint
    midpoint(s::Segment) -> MGPoint
    midpoint(coords::AbstractMatrix{T}) where {T<:Real} -> AbstractMatrix

Return the midpoint of u and v.
"""
midpoint(u::MGPoint, v::MGPoint) = (u + v) * 0.5
midpoint(s::Segment) = midpoint(s.u, s.v)
midpoint(coords::AbstractMatrix{T}
    ) where {T<:Real} = (s[1, :] + s[2, :]) * 0.5


function translate(seg, ang, dist)
    rotation = [cos(ang) -sin(ang); sin(ang) cos(ang)]
    move = rotation * toarray(unitvector(seg)) * dist
    return seg + move
end


trim_u(s::Segment, k) = Segment(s.u + (s.v - s.u) * k, s.v)
trim_v(s::Segment, k) = Segment(s.u, s.v + (s.v - s.u) * -k)
trim_uv(s::Segment, k
    ) = Segment(s.u + (s.v - s.u) * 0.5k, s.v + (s.v - s.u) * -0.5k)


cross2d(u::Point2D, v::Point2D) = u.x * v.y - u.y * v.x
cross2d(u, v) = u[1] * v[2] - u[2] * v[1]

interiorangle(u, v) = acos(dot(u, v) / (norm(u) * norm(v)))


"""
    isclockwise(vertices::AbstractMatrix{Float64}) -> Union{Bool,Nothing}

Return true/false if given vertices of a polygon in 2D space are placed clockwise/anticlockwise. Return nothing if the polygon is self-intersecting or some vertices are overlapped.  
"""
function isclockwise(vertices::AbstractMatrix{Float64})
    vlen = size(vertices, 1)
    clockwise = 0
    counter = 0
    cycpath = vcat(vertices, vertices)
    for i in 1:vlen
        p = cycpath[i, :]
        q = cycpath[i+1, :]
        r = cycpath[i+2, :]
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
        return true
    elseif round(counter / pi) == vlen - 2
        return false
    end
    return # overlapped or self-intersecting
end


function transformmatrix(scale::Point2D, rotate::Point2D, transl::Point2D)
    s = Matrix{Float64}([scale.x 0 0; 0 scale.y 0; 0 0 1])
    r = Matrix{Float64}([rotate.x -rotate.y 0; rotate.y rotate.x 0; 0 0 1])
    t = Matrix{Float64}([1 0 transl.x; 0 1 transl.y; 0 0 1])
    return t * r * s
end


function rotation(axis::Point3D, angle)
    (x1, y1, z1) = (axis.x, axis.y, axis.z)
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
