#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    Coordinates,
    cartesian2d, cartesian3d, internalcoords,
    point, segment,
    x, y, z, u, v, ux, uy, uz, vx, vy, vz,
    vector, distance,
    setcoord!,
    x_components, y_components, z_components,
    radiantophase


abstract type Coordinates end

# TODO: Matrix shape (node-wise or coordinate-wise?)


radiantophase(angle) = mod((angle + 2pi) / 2pi)
