#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export rotatable!


function rotatable!(mol; recalculate=false)
    if haskey(mol.v, :Rotatable) && !recalculate
        return
    end
    topology!(mol, recalculate=recalculate)
    elemental!(mol, recalculate=recalculate)
    # TODO: use RingBond
    pred = (i, b) -> (
        mol.v[:BondOrder][i] == 1
        && mol.v[:Degree][b.u] != 1
        && mol.v[:Degree][b.v] != 1
        && isempty(intersect(mol.v[:RingMem][b.u], mol.v[:RingMem][b.v]))
    )
    mol.v[:Rotatable] = pred.(collect(1:bondcount(mol)), mol.graph.edges)
    return
end
