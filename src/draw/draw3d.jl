#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

"""
    atom_radius(mol::SimpleMolGraph; mapping=ATOM_VANDERWAALS_RADII) -> Vector{Float64}

Return atom radii of the molecule for 3D drawing.

`mapping` is a number (same radius for all atoms) or a mapping of atomic number to radius
such as `ATOM_VANDERWAALS_RADII` and `ATOM_COVALENT_RADII`. If an element has multiple
values (`Dict{String,Float32}`), the radius of carbon is chosen by [`hybridization`](@ref)
(`"Csp3"`, `"Csp2"` or `"Csp"`) and the smallest one is chosen for the other elements.
"""
function atom_radius(mol::SimpleMolGraph; mapping=ATOM_VANDERWAALS_RADII)
    isa(mapping, Real) && return fill(Float64(mapping), nv(mol))
    desc = zeros(Float64, nv(mol))
    hyb = Symbol[]
    for i in vertices(mol)
        an = atom_number(mol[i])
        r = mapping[an]
        if isa(r, Real)
            desc[i] = r
        elseif an == 6
            # Carbon has values for each hybridization (ATOM_COVALENT_RADII)
            isempty(hyb) && (hyb = hybridization(mol))
            key = hyb[i] === :sp2 ? "Csp2" : hyb[i] === :sp ? "Csp" : "Csp3"
            desc[i] = r[key]
        else
            # For metals with low/high spin values, it's safest to choose the smallest radius
            desc[i] = minimum(values(r))
        end
    end
    return desc
end
