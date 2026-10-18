
function translate(seg::Line, ang, dist)
    # rotate the unit vector by `ang` (written out instead of a matrix product for `juliac --trim`)
    u = normalize(seg[2] - seg[1])
    move = Point2d(cos(ang) * u[1] + (-sin(ang)) * u[2], sin(ang) * u[1] + cos(ang) * u[2]) * dist
    return Line(seg[1] + move, seg[2] + move)
end


trim_u(s::Line, k) = Line(Point2d(s[1] + normalize(s[2] - s[1]) * k), Point2d(s[2]))
trim_v(s::Line, k) = Line(Point2d(s[1]), Point2d(s[2] + normalize(s[2] - s[1]) * -k))
trim_uv(s::Line, k) = Line(Point2d(s[1] + normalize(s[2] - s[1]) * k), Point2d(s[2] + normalize(s[2] - s[1]) * -k))

# TODO: check cases abs(dot(u, v) / (norm(u) * norm(v))) > 1
interiorangle(u, v) = acos(max(-1, min(1, dot(u, v) / (norm(u) * norm(v)))))


"""
    isclockwise(vertices::Vector{Point2d}) -> Union{Bool,Nothing}

Return true/false if given vertices of a polygon in 2D space are placed clockwise/anticlockwise. Return nothing if the polygon is self-intersecting or some vertices are overlapped.  
"""
function isclockwise(vertices::Vector{Point2d})
    vlen = length(vertices)
    clockwise = 0
    counter = 0
    cycpath = vcat(vertices, vertices)
    for i in 1:vlen
        p = cycpath[i]
        q = cycpath[i + 1]
        r = cycpath[i + 2]
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
    return # overlapped or self-intersecting
end


# 3x3 affine matrix [a11 a12 a13; a21 a22 a23; 0 0 1]
# (matrix literals are built by Vararg hvcat, which is not compilable with `juliac --trim`)
function affine3(a11, a12, a13, a21, a22, a23)
    m = zeros(Float64, 3, 3)
    m[1, 1], m[1, 2], m[1, 3] = a11, a12, a13
    m[2, 1], m[2, 2], m[2, 3] = a21, a22, a23
    m[3, 3] = 1.0
    return m
end

function transformmatrix(scale::Point2d, rotate::Point2d, transl::Point2d)
    s = affine3(scale[1], 0.0, 0.0, 0.0, scale[2], 0.0)
    r = affine3(rotate[1], -rotate[2], 0.0, rotate[2], rotate[1], 0.0)
    t = affine3(1.0, 0.0, transl[1], 0.0, 1.0, transl[2])
    return t * r * s
end


function transformmatrix(seg::Line, xscale::Float64, yscale::Float64)
    dist = norm(seg[2] - seg[1])
    scalef = Point2d(dist * xscale, yscale)
    rotatef = normalize(seg[2] - seg[1])
    translf = seg[1]
    return transformmatrix(scalef, rotatef, translf)
end