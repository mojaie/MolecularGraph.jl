#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    mcssize,
    mcsview



"""
    mcssize(mol1::GraphMol, mol2::GraphMol, subgraphtype=:EdgeInduced;
            atommatcher=atommatch, bondmatcher=bondmatch,
            timeout=60, kwargs...) -> Int

Compute maximum common substructure (MCS) and return number of nodes
(subgraphtype=:NodeInduced) or edges (subgraphtype=:EdgeInduced).

# Options

- atommatcher(Function): customized `atommatcher` function
- bondmatcher(Function): customized `bondmatcher` function
- timeout(Int): abort calculation (default 60 seconds)
- algorithm(Symbol): `:Clique`
- c_clique_constraint(Bool): if true, calculate c-clique

"""
function mcssize(mol1, mol2; subgraphtype=:EdgeInduced,
                 atommatcher=atommatch, bondmatcher=bondmatch,
                 timeout=60, kwargs...)
    afunc = atommatcher(mol1, mol2)
    bfunc = bondmatcher(mol1, mol2)
    mcs = maximumcommonsubgraph(
        mol1, mol2, subgraphtype=subgraphtype,
        nodematcher=afunc, edgematcher=bfunc, timeout=timeout; kwargs...)
    return length(mcs)
end


"""
    mcsmol(mol1::GraphMol, mol2::GraphMol; mode=:MCES) -> MapMol

Return maximum common substructure (MCS) as MolGraphView.
"""
function mcsview(mol1, mol2; kwargs...)
    mcs = maximumcommonsubgraph(mol1, mol2; kwargs...)
end
