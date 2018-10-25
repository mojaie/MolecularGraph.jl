#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    vec2d,
    vec3d,
    posX,
    posY,
    posZ,
    interiorangle,
    vecpair,
    utov,
    swap,
    vecU,
    vecV,
    translate,
    trimU,
    trimV,
    trimUV,
    rotation,
    isclockwise,
    radiantophase,
    transformmatrix

import LinearAlgebra: cross


# 2D vector operation

vec2d(x::AbstractFloat, y::AbstractFloat) = SVector{2}(x, y)
vec2d(v::AbstractArray) = vec2d(v[1], v[2])
posX(v::SVector{2}) = v[1]
posY(v::SVector{2}) = v[2]
cross(u::SVector{2}, v::SVector{2}) = u[1] * v[2] - u[2] * v[1]
interiorangle(u::SVector{2}, v::SVector{2}) = acos(
    dot(u, v) / (norm(u) * norm(v)))
vecpair(u::SVector{2}, v::SVector{2}) = SMatrix{2,2}(u[1], v[1], u[2], v[2])



# 2D vector pair operation

utov(uv::SMatrix{2,2}) = uv[2, :] - uv[1, :]
swap(uv::SMatrix{2,2}) = SMatrix{2,2}(uv[2], uv[1], uv[4], uv[3])
vecU(uv::SMatrix{2,2}) = SVector{2}(uv[1, :])
vecV(uv::SMatrix{2,2}) = SVector{2}(uv[2, :])


function translate(uv::SMatrix{2,2}, angle::Real, dist::Real)
    rot = rotation(angle)
    move = rot * normalize(utov(uv)) * dist
    uv .+ transpose(move)
end


trimU(uv::SMatrix{2,2}, k::Real) = uv + [transpose(utov(uv) * k); 0 0]
trimV(uv::SMatrix{2,2}, k::Real) = uv + [0 0; transpose(utov(uv) * -k)]
trimUV(uv::SMatrix{2,2}, k::Real) = uv + [
    transpose(utov(uv) * k / 2); transpose(utov(uv) * -k / 2)]



# 3D vector operation

vec3d(x::AbstractFloat, y::AbstractFloat,
      z::AbstractFloat) = SVector{3}(x, y, z)
vec3d(v::AbstractArray) = vec3d(v[1], v[2], v[3])
posX(v::SVector{3}) = v[1]
posY(v::SVector{3}) = v[2]
posZ(v::SVector{3}) = v[3]



# Rotation matrix

rotation(a::Real) = SMatrix{2,2}([cos(a) -sin(a); sin(a) cos(a)])

function rotation(axis::SVector{3}, angle::Real)
    (x, y, z) = axis
    c = cos(angle)
    s = sin(angle)
    a12 = x * y * (1 - c)
    a13 = x * z * (1 - c)
    a23 = y * z * (1 - c)
    SMatrix{3,3}([
        (c + x^2 * (1 - c)) (a12 - z * s) (a13 + y * s);
        (a12 + z * s) (c + y^2 * (1 - c)) (a23 - x * s);
        (a13 - y * s) (a23 + x * s) (c + z^2 * (1 - c))
    ])
end



# Convenient methods

function isclockwise(vertices::AbstractArray{SVector{2}})
    vlen = length(vertices)
    clockwise = 0
    counter = 0
    cycver = cat(vertices, vertices, dims=1)
    for i = 1:vlen
        (p, q, r) = cycver[i:i+2]
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
        true
    elseif round(counter / pi) == vlen - 2
        false
    else
        NaN
    end
end


function radiantophase(angle::Real)
    mod((angle + 2pi) / 2pi)
end


function transformmatrix(scale::SVector{2}, rotvec::SVector{2},
                         transl::SVector{2})
    scale = SMatrix{3,3}([scale[1] 0 0; 0 scale[2] 0; 0 0 1])
    rot = SMatrix{3,3}([rotvec[1] -rotvec[2] 0; rotvec[2] rotvec[1] 0; 0 0 1])
    tl = SMatrix{3,3}([1 0 transl[1]; 0 1 transl[2]; 0 0 1])
    tl * rot * scale
end
