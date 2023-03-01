#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    coords2d, sdfcoords2d, sdfcoords2d!, coordgen, coordgen!

using coordgenlibs_jll


coords2d(mol::MolGraph) = get_descriptor(mol, :v_coords2d)


function sdfcoords2d(mol::MolGraph)
    coords = zeros(Float64, nv(mol), 2)
    for i in vertices(mol)
        coords[i, :] = get_prop(mol, i, :coords)[1:2]
    end
    return coords
end
sdfcoords2d!(mol::MolGraph) = set_descriptor!(mol, :v_coords2d, sdfcoords2d(mol))




"""
    coordgen(mol::MolGraph) -> Tuple{Array{Int,1},Array{Int,1}}

Generate 2D coordinates by using Schrodinger's coordgenlibs.

This will returns a tuple of `coords` and `styles` arrays. `coords` is a size(n, 2) matrix where n is atom count, which stores 2D coordinates (x, y) of each atoms. `styles` is a size e vector of wedge notation of stereobond, where e is bond count.
"""
function coordgen(mol::MolGraph)
    minmol = ccall((:getSketcherMinimizer, libcoordgen), Ptr{Cvoid}, ())
    atoms = Ptr{Cvoid}[]
    bonds = Ptr{Cvoid}[]

    # Atoms
    for a in atom_symbol(mol)
        atom = ccall(
            (:setAtom, libcoordgen), Ptr{Cvoid},
            (Ptr{Cvoid}, Int), minmol, atomnumber(a)
        )
        push!(atoms, atom)
    end

    # Bonds
    bondorder_ = bond_order(mol)
    for (i, e) in enumerate(edges(mol))
        bond = ccall(
            (:setBond, libcoordgen), Ptr{Cvoid},
            (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Int),
            minmol, atoms[src(e)], atoms[dst(e)], bondorder_[i]
        )
        push!(bonds, bond)
    end

    ccall(
        (:assignBondsAndNeighbors, libcoordgen), Cvoid,
        (Ptr{Cvoid},), minmol
    )

    # Stereocenter
    if has_prop(mol, :stereocenter)
        for (n, stereo) in get_prop(mol, :stereocenter)
            ccall(
                (:setStereoCenter, libcoordgen), Cvoid,
                (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Int),
                atoms[n], atoms[stereo[1]], atoms[stereo[2]], atoms[stereo[3]], stereo[4]
            )
        end
    end

    # Stereobond
    if has_prop(mol, :stereobond)
        for (e, stereo) in get_prop(mol, :stereobond)
            ccall(
                (:setStereoBond, libcoordgen), Cvoid,
                (Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}, Int),
                bonds[edge_rank(mol, e)], atoms[stereo[1]], atoms[stereo[2]], stereo[3]
            )
        end
    end

    # Optimize
    ccall(
        (:runGenerateCoordinates, libcoordgen), Cvoid,
        (Ptr{Cvoid},), minmol
    )

    # Output
    coords = zeros(Float64, nv(mol), 2)
    for i in vertices(mol)
        px = ccall(
            (:getAtomX, libcoordgen), Float32, (Ptr{Cvoid},), atoms[i])
        py = ccall(
            (:getAtomY, libcoordgen), Float32, (Ptr{Cvoid},), atoms[i])
        coords[i, :] = [px, py]
    end
    styles = Int[]
    for e in edges(mol)
        hasstereo = ccall(
            (:hasStereochemistryDisplay, libcoordgen), Bool,
            (Ptr{Cvoid},), bonds[edge_rank(mol, e)])
        if !hasstereo
            push!(styles, 0)
            continue
        end
        iswedge = ccall(
            (:isWedge, libcoordgen), Bool, (Ptr{Cvoid},), bonds[edge_rank(mol, e)])
        isrev = ccall(
            (:isReversed, libcoordgen), Bool, (Ptr{Cvoid},), bonds[edge_rank(mol, e)])
        push!(styles, (iswedge ? 1 : 6) + (isrev ? 1 : 0))
    end
    return coords, styles
end

function coordgen!(mol::MolGraph)
    coords, styles = coordgen(mol)
    set_descriptor!(mol, :v_coords2d, coords)
    set_descriptor!(mol, :v_styles, styles)
    return
end