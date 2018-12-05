#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    default_annotation!


function default_annotation!(mol::VectorMol)
    molgraph_topology!(mol)
    default_annot!(mol)
    return
end


# Graph based stereoisomer determination
#:BondEorZ
#    -> SMILES_CisTrans
#:CIP_Rule  # TODO Dont use CIP for chirality info. use Index like SMILES
#    -> SDFile_BondType
#    -> SMILES_Stereo
