#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#


abstract type MGPoint end  # TODO: to be removed, replaced by GeometryBasics
abstract type Coordinates end


radiantophase(angle) = mod((angle + 2pi) / 2pi)
