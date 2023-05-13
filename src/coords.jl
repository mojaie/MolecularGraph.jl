#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    has_coords, sdf_coords2d, coords2d, coords3d, coordgen, coordgen!

using coordgenlibs_jll


function has_coords(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_coords2d) && return true
    hasfield(vproptype(mol), :coords) || return false
    for i in vertices(mol)
        isnothing(get_prop(mol, i, :coords)) && return false
    end
    return true
end

function sdf_coords2d(mol::SimpleMolGraph)
    coords = zeros(Float64, nv(mol), 2)
    for i in vertices(mol)
        coords[i, :] = get_prop(mol, i, :coords)[1:2]
    end
    return coords
end

function coords2d(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_coords2d) && return get_cache(mol, :v_coords2d)
    has_coords(mol) || error("no coordinates. use coordgen")
    return sdf_coords2d(mol)
end


function coords3d(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    has_cache(mol, :v_coords3d) && return get_cache(mol, :v_coords3d)
    has_coords(mol) || error("no coordinates. use coordgen")
    coords = zeros(Float64, nv(mol), 3)
    for i in vertices(mol)
        coords[i, :] = get_prop(mol, i, :coords)[1:3]
    end
    return coords
end


"""
    coordgen(mol::MolGraph) -> Tuple{Array{Int,1},Array{Int,1}}

Generate 2D coordinates by using Schrodinger's coordgenlibs.

This will returns a tuple of `coords` and `styles` arrays. `coords` is a size(n, 2) matrix
where n is atom count, which stores 2D coordinates (x, y) of each atoms.
`styles` is a size e vector of wedge notation of stereobond, where e is bond count.
"""
function coordgen(g, atomsymbol_, bondorder_, stereocenters, stereobonds)
    minmol = @ccall libcoordgen.getSketcherMinimizer()::Ptr{Cvoid}
    atoms = Ptr{Cvoid}[]
    bonds = Ptr{Cvoid}[]
    er = Dict(e => i for (i, e) in enumerate(edges(g)))

    # Atoms
    for a in atomsymbol_
        atom = @ccall libcoordgen.setAtom(
            minmol::Ptr{Cvoid}, atomnumber(a)::Cint)::Ptr{Cvoid}
        push!(atoms, atom)
    end

    # Bonds
    for (i, e) in enumerate(edges(g))
        bond = @ccall libcoordgen.setBond(
            minmol::Ptr{Cvoid}, atoms[src(e)]::Ptr{Cvoid},
            atoms[dst(e)]::Ptr{Cvoid}, bondorder_[i]::Cint)::Ptr{Cvoid}
        push!(bonds, bond)
    end

    @ccall libcoordgen.assignBondsAndNeighbors(minmol::Ptr{Cvoid})::Cvoid

    # Stereocenter
    for (n, stereo) in stereocenters
        @ccall libcoordgen.setStereoCenter(
            atoms[n]::Ptr{Cvoid}, atoms[stereo[1]]::Ptr{Cvoid},
            atoms[stereo[2]]::Ptr{Cvoid}, atoms[stereo[3]]::Ptr{Cvoid}, stereo[4]::Cint)::Cvoid
    end

    # Stereobond
    for (e, stereo) in stereobonds
        @ccall libcoordgen.setStereoBond(
            bonds[er[e]]::Ptr{Cvoid}, atoms[stereo[1]]::Ptr{Cvoid},
            atoms[stereo[2]]::Ptr{Cvoid}, stereo[3]::Cint
        )::Cvoid
    end

    # Optimize
    @ccall libcoordgen.runGenerateCoordinates(minmol::Ptr{Cvoid})::Cvoid

    # Output
    coords = zeros(Float64, nv(g), 2)
    for i in vertices(g)
        px = @ccall libcoordgen.getAtomX(atoms[i]::Ptr{Cvoid})::Cfloat
        py = @ccall libcoordgen.getAtomY(atoms[i]::Ptr{Cvoid})::Cfloat
        coords[i, :] = [px, py]
    end
    styles = Vector{Symbol}(undef, ne(g))
    for (i, e) in enumerate(edges(g))
        hasstereo = @ccall libcoordgen.hasStereochemistryDisplay(
            bonds[i]::Ptr{Cvoid})::Bool
        if !hasstereo
            styles[i] = :none
            continue
        end
        iswedge = @ccall libcoordgen.isWedge(bonds[i]::Ptr{Cvoid})::Bool
        isrev = @ccall libcoordgen.isReversed(bonds[i]::Ptr{Cvoid})::Bool
        if iswedge
            styles[i] = isrev ? :revup : :up
        else
            styles[i] = isrev ? :revdown : :down
        end
    end
    # TODO: keep wave bond in SDFile
    return coords, styles
end


function coordgen(mol::SimpleMolGraph)
    get_state(mol, :has_updates) && dispatch!(mol, :updater)
    return coordgen(
        mol.graph, atom_symbol(mol), bond_order(mol),
        get_prop(mol, :stereocenter), get_prop(mol, :stereobond)
    )
end


function coordgen!(mol::SimpleMolGraph)
    coords, styles = coordgen(
        mol.graph, atom_symbol(mol), bond_order(mol),
        get_prop(mol, :stereocenter), get_prop(mol, :stereobond)
    )
    set_cache!(mol, :v_coords2d, coords)
    set_cache!(mol, :e_coordgen_bond_style, styles)
end