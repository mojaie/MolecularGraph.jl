#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    mcismol, mcismolsize,
    mcesmol, mcesmolsize


"""
    mcismol(mol1::UndirectedGraph, mol2::UndirectedGraph; kwargs...
        ) -> Tuple{Dict{Int,Int},Symbol}

Compute maximum common induced substructure (MCIS) of mol1 and mol2.

## Keyword arguments

- connected(Bool): if true, apply connected MCS constraint.
- topological(Bool): if true, apply topological constraint.
- diameter(Int): distance cutoff for topological constraint.
- tolerance(Int): distance mismatch tolerance for topological constraint.
- timeout(Int): abort calculation and return suboptimal results if the execution
time has reached the given value (default=60, in seconds).
- targetsize(Int): abort calculation and return suboptimal result so far if the
given mcs size achieved.

# References

1. Kawabata, T. (2011). Build-Up Algorithm for Atomic Correspondence between
Chemical Structures. Journal of Chemical Information and Modeling, 51(8),
1775–1787. https://doi.org/10.1021/ci2001023
1. https://www.jstage.jst.go.jp/article/ciqs/2017/0/2017_P4/_article/-char/en
"""
function mcismol(mol1::UndirectedGraph, mol2::UndirectedGraph; kwargs...)
    afunc = atommatch(mol1, mol2)
    bfunc = bondmatch(mol1, mol2)
    return findmcis(
        mol1, mol2, nodematcher=afunc, edgematcher=bfunc; kwargs...)
end

mcismolsize(mol1, mol2; kwargs...) = length(mcismol(mol1, mol2; kwargs...)[1])



"""
    mcesmol(mol1::UndirectedGraph, mol2::UndirectedGraph; kwargs...
        ) -> Tuple{Dict{Int,Int},Symbol}

Compute maximum common edge induced substructure (MCES) of mol1 and mol2.
"""
function mcesmol(mol1::UndirectedGraph, mol2::UndirectedGraph; kwargs...)
    afunc = atommatch(mol1, mol2)
    bfunc = bondmatch(mol1, mol2)
    return findmces(
        mol1, mol2, nodematcher=afunc, edgematcher=bfunc; kwargs...)
end

mcesmolsize(mol1, mol2; kwargs...) = length(mcesmol(mol1, mol2; kwargs...)[1])
