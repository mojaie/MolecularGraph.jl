#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    Coordinates, Cartesian2D, Cartesian3D, InternalCoords


abstract type Coordinates end

# TODO: Matrix shape (node-wise or coordinate-wise?)

mutable struct Cartesian2D <: Coordinates
    coords::Matrix{Float64}
end


mutable struct Cartesian3D <: Coordinates
    coords::Matrix{Float64}
end


mutable struct InternalCoords <: Coordinates
    labels::Matrix{Union{Int,Nothing}} # atom1, atom2 and atom3
    geometry::Matrix{Union{Float64,Nothing}} # distance, angle and dihedral
    nodekeys::Dict{Int,Int} # labels, graph node keys
end
