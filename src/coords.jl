#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

to_dict(
    ::Val{:default}, ::Val{:coords2d}, gprop::MolProperty
) = [[collect(cd) for cd in cds] for cds in gprop.coords2d]
reconstruct(
    ::Val{:coords2d}, gprop::MolProperty, data
) = [[Point2d(cd...) for cd in cds] for cds in data]
remap(  # TODO
    ::Val{:coords2d}, gprop::MolProperty{T}, vmap::Dict{T,T}
) where T = gprop.coords2d

to_dict(
    ::Val{:default}, ::Val{:coords3d}, gprop::MolProperty
) = [[collect(cd) for cd in cds] for cds in gprop.coords3d]
reconstruct(
    ::Val{:coords3d}, gprop::MolProperty, data
) = [[Point3d(cd...) for cd in cds] for cds in data]
remap(  # TODO
    ::Val{:coords3d}, gprop::MolProperty{T}, vmap::Dict{T,T}
) where T = gprop.coords3d

to_dict(
    ::Val{:default}, ::Val{:draw2d_bond_style}, gprop::MolProperty
) = [[string(s) for s in sty] for sty in gprop.draw2d_bond_style]
reconstruct(
    ::Val{:draw2d_bond_style}, gprop::MolProperty, data
) = [[Symbol(s) for s in sty] for sty in data]
remap(  # TODO
    ::Val{:draw2d_bond_style}, gprop::MolProperty{T}, vmap::Dict{T,T}
) where T = gprop.draw2d_bond_style


function coords_from_sdf!(mol)
    # Initializer. should be called before stereochem initializers.
    zrange = nv(mol) == 0 ? (0, 0) : extrema(get_prop(mol, i, :coords)[3] for i in vertices(mol))
    # Embed 3D coords to 2D
    push!(mol.gprops.coords2d, [Point2d(get_prop(mol, i, :coords)[1:2]...) for i in vertices(mol)])
    if zrange[2] - zrange[1] > 0.001  # 3D coords available
        push!(mol.gprops.coords3d, [Point3d(get_prop(mol, i, :coords)[1:3]...) for i in vertices(mol)])
    end
    # Bond style in 2D notation (wedges)
    bondorder = [get_prop(mol, e, :order) for e in edges(mol)]
    bondnotation = [get_prop(mol, e, :notation) for e in edges(mol)]
    isordered = [get_prop(mol, e, :isordered) for e in edges(mol)]
    arr = Vector{Symbol}(undef, length(bondorder))
    for i in 1:length(bondorder)
        if bondnotation[i] == 3
            arr[i] = :cis_trans
        elseif bondorder[i] != 1
            arr[i] = :none
        elseif bondnotation[i] == 1
            arr[i] = isordered[i] ? :up : :revup
        elseif bondnotation[i] == 6
            arr[i] = isordered[i] ? :down : :revdown
        elseif bondnotation[i] == 4
            arr[i] = :unspecified
        else
            arr[i] = :none
        end
    end
    push!(mol.gprops.draw2d_bond_style, arr)
end


function coords2d(mol::SimpleMolGraph, i::Int)
    dispatch_update!(mol)
    return mol.gprops.coords2d[i]
end
coords2d(mol) = coords2d(mol, 1)
has_coords2d(mol) = length(mol.gprops.coords2d) > 0


function coords3d(mol::SimpleMolGraph, i::Int)
    dispatch_update!(mol)
    return mol.gprops.coords3d[i]
end
coords3d(mol) = coords3d(mol, 1)
has_coords3d(mol) = length(mol.gprops.coords3d) > 0


function draw2d_bond_style(mol::SimpleMolGraph, i::Int)
    dispatch_update!(mol)
    return mol.gprops.draw2d_bond_style[i]
end
draw2d_bond_style(mol) = draw2d_bond_style(mol, 1)


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
            minmol::Ptr{Cvoid}, atom_number(a)::Cint)::Ptr{Cvoid}
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
    coords = Vector{Point2d}(undef, nv(g))
    for i in vertices(g)
        px = @ccall libcoordgen.getAtomX(atoms[i]::Ptr{Cvoid})::Cfloat
        py = @ccall libcoordgen.getAtomY(atoms[i]::Ptr{Cvoid})::Cfloat
        coords[i] = Point2d(px, py)
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
    dispatch_update!(mol)
    return coordgen(
        mol.graph, atom_symbol(mol), bond_order(mol),
        mol.gprops.stereocenter, mol.gprops.stereobond
    )
end


function coordgen!(mol::SimpleMolGraph)
    # Initialize
    empty!(mol.gprops.coords2d)
    empty!(mol.gprops.draw2d_bond_style)
    # TODO: unspecified stereochem in SMILES
    coords, styles = coordgen(
        mol.graph, atom_symbol(mol), bond_order(mol),
        get_prop(mol, :stereocenter), get_prop(mol, :stereobond)
    )
    push!(mol.gprops.coords2d, coords)
    push!(mol.gprops.draw2d_bond_style, styles)
end