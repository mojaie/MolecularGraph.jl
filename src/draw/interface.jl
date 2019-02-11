#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    singlebond!,
    wedged!,
    dashedwedged!,
    wavesingle!,
    doublebond!,
    clockwisedouble!,
    counterdouble!,
    crossdouble!,
    triplebond!,
    atomsymbolright!,
    atomsymbolcenter!,
    atomsymbolleft!,
    atom_annotation!


using MolecularGraph.Geometry: _point, _vector, _u, _v


# Required functions for drawing canvas
# Single bonds
function singlebond! end
function wedged! end
function dashedwedged! end
function wavesingle! end

# Multiple bonds
function doublebond! end
function clockwisedouble! end
function counterdouble! end
function crossdouble! end
function triplebond! end

# Atom symbols
function atomsymbolright! end
function atomsymbolcenter! end
function atomsymbolleft! end
function atom_annotation! end
