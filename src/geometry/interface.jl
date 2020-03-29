#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    Point, Coordinates,
    radiantophase


abstract type Point end
abstract type Coordinates end


radiantophase(angle) = mod((angle + 2pi) / 2pi)
