#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    internalcoords,
    label1, label2, label3, distance, angle, dihedral,
    setcoords!


# TODO: Gaussian compatible z-matrix format


function internalcoords(labels, geometry)
    size(labels, 2) == 3 || throw(
        DimensionMismatch("Unexpected matrix size $(size(coords))"))
    size(geometry, 2) == 3 || throw(
        DimensionMismatch("Unexpected matrix size $(size(coords))"))
    return InternalCoords(labels, geometry, Dict())
end

function internalcoords(size::Int)
    labels = fill(nothing, (size, 3))
    geometry = fill(nothing, (size, 3))
    return InternalCoords(labels, geometry, Dict())
end


Base.length(coords::InternalCoords) = size(coords.labels, 1)


label1(coords::InternalCoords, i::Int) = coords.labels[i, 1]
label2(coords::InternalCoords, i::Int) = coords.labels[i, 2]
label3(coords::InternalCoords, i::Int) = coords.labels[i, 3]
distance(coords::InternalCoords, i::Int) = coords.geometry[i, 1]

_angle(coords::InternalCoords, i::Int) = coords.geometry[i, 2]
Base.angle(coords::InternalCoords, i::Int) = _angle(coords, i) * pi

_dihedral(coords::InternalCoords, i::Int) = coords.geometry[i, 3]
dihedral(coords::InternalCoords, i::Int) = _dihedral(coords, i) * pi


function setcoord!(coords::InternalCoords, i::Int, labels, geometry)
    coords.labels[i, :] = labels
    coords.geometry[i, :] = geometry
end
