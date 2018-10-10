#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#


function draw(canvas::Canvas, mol::MolecularGraph)
    required_descriptor("ScaleAndCenter")
    mlb = mol.size2d[3]
    if length(atomvector(mol)) == 0
        return
    end
    
end
