#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#


export
    molweight,
    H_acceptor_count,
    H_donor_count,
    wclogp,
    rotatable_count


"""
    molweight(mol::VectorMol; digits=2) -> Float64

Return standard molecular weight.
"""
function molweight(mol::VectorMol; digits=2)
    return round(reduce(+, mol.v[:MolWeight]; init=0), digits=digits)
end


"""
    H_acceptor_count(mol::VectorMol) -> Int

Return the number of hydrogen bond acceptors (N, O and F).
"""
function H_acceptor_count(mol::VectorMol)
    return reduce(+, mol.v[:H_Acceptor]; init=0)
end


"""
    H_donor_count(mol::VectorMol) -> Int

Return the number of hydrogen bond donors (O and N attached to hydrogens).
"""
function H_donor_count(mol::VectorMol)
    return reduce(+, mol.v[:H_Donor]; init=0)
end


"""
    wclogp(mol::VectorMol) -> Float64

Return predicted logP value calculated by using Wildman and Crippen method.

# Reference

1. Wildman, S. A. and Crippen, G. M. (1999). Prediction of Physicochemical
Parameters by Atomic Contributions. Journal of Chemical Information and
Modeling, 39(5), 868–873. https://doi.org/10.1021/ci990307l
"""
function wclogp(mol::VectorMol; digits=2)
    return round(reduce(+, mol.v[:WCLogPContrib]; init=0), digits=digits)
end


"""
    rotatable_count(mol::VectorMol) -> Int

Return the number of rotatable bonds.
"""
function rotatable_count(mol::VectorMol)
    return reduce(+, mol.v[:Rotatable]; init=0)
end
