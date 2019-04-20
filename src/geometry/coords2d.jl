#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    Cartesian2D, Point2D, Segment2D,
    interiorangle, translate,
    midpoint, trim_u, trim_v, trim_uv,
    transformmatrix, isclockwise


struct Cartesian2D{T<:AbstractMatrix{Float64}} <: Coordinates
    coords::T
end


struct Segment2D{T<:AbstractMatrix{Float64}}
    coords::T
end


struct Point2D{T<:AbstractVector{Float64}}
    pos::T
end


cartesian2d(data::AbstractMatrix{Float64}) = Cartesian2D(data)
cartesian2d(c2d::Cartesian2D, vertices::Vector{Int}
    ) = Cartesian2D(@view c2d.coords[vertices, :])

segment(c2d::Cartesian2D, u, v) = Segment2D(@view c2d.coords[[u, v], :])
segment(u::Point2D, v::Point2D) = Segment2D(transpose(hcat(u.pos, v.pos)))

point(c2d::Cartesian2D, i::Int) = Point2D(@view c2d.coords[i, :])
point(x::Real, y::Real) = Point2D([float(x), float(y)])

Base.length(c2d::Cartesian2D) = size(c2d.coords, 1)

Base.:+(u::Point2D, v::Point2D) = Point2D(u.pos + v.pos)
Base.:-(u::Point2D, v::Point2D) = Point2D(u.pos - v.pos)
Base.:+(u::Point2D, v::Tuple{Real,Real}) = Point2D(u.pos + collect(v))
Base.:-(u::Point2D, v::Tuple{Real,Real}) = Point2D(u.pos - collect(v))
Base.:*(p::Point2D, f::Float64) = Point2D(p.pos * f)

Base.:+(seg::Segment2D, t::Point2D) = Segment2D(seg.coords .+ transpose(t.pos))
Base.:-(seg::Segment2D, t::Point2D) = Segment2D(seg.coords .- transpose(t.pos))
Base.:*(seg::Segment2D, f::Float64) = Segment2D(seg.coords * f)

LinearAlgebra.norm(p::Point2D) = norm(p.pos)
LinearAlgebra.normalize(p::Point2D) = Point2D(normalize(p.pos))
LinearAlgebra.dot(u::Point2D, v::Point2D) = dot(u.pos, v.pos)
LinearAlgebra.cross(u::Point2D, v::Point2D) = x(u) * y(v) - y(u) * x(v)


x(point::Point2D) = point.pos[1]
y(point::Point2D) = point.pos[2]

u(segment::Segment2D) = Point2D(@view segment.coords[1, :])
v(segment::Segment2D) = Point2D(@view segment.coords[2, :])
ux(segment::Segment2D) = segment.coords[1]
uy(segment::Segment2D) = segment.coords[3]
vx(segment::Segment2D) = segment.coords[2]
vy(segment::Segment2D) = segment.coords[4]



vector(segment::Segment2D) = v(segment) - u(segment)
distance(segment::Segment2D) = norm(vector(segment))
interiorangle(u::Point2D, v::Point2D) = acos(dot(u, v) / (norm(u) * norm(v)))


function translate(s::Segment2D, ang, dist)
    rotation = [cos(ang) -sin(ang); sin(ang) cos(ang)]
    move = Point2D(rotation * normalize(vector(s)).pos * dist)
    return s + move
end


midpoint(s::Segment2D) = (u(s) + v(s)) * 0.5

trim_u(s::Segment2D, k) = segment(u(s) + vector(s) * k, v(s))
trim_v(s::Segment2D, k) = segment(u(s), v(s) + vector(s) * -k)
trim_uv(s::Segment2D, k
    ) = segment(u(s) + vector(s) * 0.5k, v(s) + vector(s) * -0.5k)


function setcoord!(segment::Segment2D, i::Int, point::Point2D)
    segment.coords[i, :] = point.pos
    return
end


x_components(c2d::Cartesian2D) = @view c2d.coords[:, 1]
y_components(c2d::Cartesian2D) = @view c2d.coords[:, 2]



function transformmatrix(scale::Point2D, rotate::Point2D, translate::Point2D)
    s = Matrix{Float64}([x(scale) 0 0; 0 y(scale) 0; 0 0 1])
    r = Matrix{Float64}([x(rotate) -y(rotate) 0; y(rotate) x(rotate) 0; 0 0 1])
    t = Matrix{Float64}([1 0 x(translate); 0 1 y(translate); 0 0 1])
    return t * r * s
end


function isclockwise(vertices::Cartesian2D)
    vlen = length(vertices)
    clockwise = 0
    counter = 0
    cycpath = vcat(vertices.coords, vertices.coords)
    for i in 1:vlen
        p = Point2D(cycpath[i, :])
        q = Point2D(cycpath[i+1, :])
        r = Point2D(cycpath[i+2, :])
        cp = cross(q - p, r - p)
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
    return # overlapped
end
