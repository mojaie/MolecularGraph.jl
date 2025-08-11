#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

# TODO: just remap nodes if still all existing vertices have coords.

reconstruct(::Val{:coords2d}, ::Type{T}, @nospecialize(data)
    ) where T <: AbstractProperty = [[Point2d(cd...) for cd in cds] for cds in data]
to_dict(
    ::Val{:coords2d}, ::Val{:default}, gprop::AbstractProperty
) = [[collect(cd) for cd in cds] for cds in gprop.coords2d]

function remap!(
        ::Val{:coords2d}, desc::SimpleMolProperty{T}, vmap::Vector{T},
        edges::Vector{Edge{T}}) where T <: Integer
    revv = Dict(v => i for (i, v) in enumerate(vmap))
    container = Vector{Vector{Point2d}}(undef, length(desc.coords2d))
    for (i, oldcds) in enumerate(desc.coords2d)
        newcds = Vector{Point2d}(undef, length(revv))
        for (oldp, newp) in revv
            newcds[newp] = oldcds[oldp]
        end
        container[i] = newcds
    end
    empty!(desc.coords2d)
    append!(desc.coords2d, container)
    return
end

reconstruct(::Val{:coords3d}, ::Type{T}, @nospecialize(data)
    ) where T <: AbstractProperty = [[Point3d(cd...) for cd in cds] for cds in data]
to_dict(
    ::Val{:coords3d}, ::Val{:default}, gprop::AbstractProperty
) = [[collect(cd) for cd in cds] for cds in gprop.coords3d]

function remap!(
        ::Val{:coords3d}, desc::SimpleMolProperty{T}, vmap::Vector{T},
        edges::Vector{Edge{T}}) where T <: Integer
    revv = Dict(v => i for (i, v) in enumerate(vmap))
    container = Vector{Vector{Point3d}}(undef, length(desc.coords3d))
    for (i, oldcds) in enumerate(desc.coords3d)
        newcds = Vector{Point3d}(undef, length(revv))
        for (oldp, newp) in revv
            newcds[newp] = oldcds[oldp]
        end
        container[i] = newcds
    end
    empty!(desc.coords3d)
    append!(desc.coords3d, container)
    return
end

reconstruct(::Val{:draw2d_bond_style}, ::Type{T}, @nospecialize(data)
    ) where T <: AbstractProperty = [[Symbol(s) for s in sty] for sty in data]
to_dict(
    ::Val{:draw2d_bond_style}, ::Val{:default}, gprop::AbstractProperty
) = [[string(s) for s in sty] for sty in gprop.draw2d_bond_style]

function remap!(
        ::Val{:draw2d_bond_style}, desc::SimpleMolProperty{T}, vmap::Vector{T},
        edges::Vector{Edge{T}}) where T <: Integer
    container = Vector{Vector{Symbol}}(undef, length(desc.draw2d_bond_style))
    revv = Dict(v => i for (i, v) in enumerate(vmap))
    newedges = sort([u_edge(T, revv[src(e)], revv[dst(e)]) for e in edges
            if src(e) in vmap && dst(e) in vmap])
    olderank = Dict(e => i for (i, e) in enumerate(edges))
    emap = Dict(i => olderank[u_edge(T, vmap[src(e)], vmap[dst(e)])] for (i, e) in enumerate(newedges))
    inv = Dict{Symbol,Symbol}(:up => :revup, :revup => :up, :down => :revdown, :revdown => :down)
    for (i, styles) in enumerate(desc.draw2d_bond_style)
        cont = Vector{Symbol}(undef, length(emap))
        for (j, e) in enumerate(newedges)
            dr = styles[emap[j]]
            if dr in [:up, :revup, :down, :revdown] && (src(e) < dst(e)) !== (vmap[src(e)] < vmap[dst(e)])
                cont[j] = inv[dr]
            else
                cont[j] = dr
            end
        end
        container[i] = cont
    end
    empty!(desc.draw2d_bond_style)
    append!(desc.draw2d_bond_style, container)
    return
end


function coords_from_sdf!(mol::SimpleMolGraph)
    # Initializer. should be called before stereochem initializers.
    zrange = nv(mol) == 0 ? (0, 0) : extrema(get_prop(mol, i, :coords)[3] for i in vertices(mol))
    # Embed 3D coords to 2D
    push!(
        mol.gprops.descriptors.coords2d,
        [Point2d(get_prop(mol, i, :coords)[1:2]...) for i in vertices(mol)])
    if zrange[2] - zrange[1] > 0.001  # 3D coords available
        push!(
            mol.gprops.descriptors.coords3d,
            [Point3d(get_prop(mol, i, :coords)[1:3]...) for i in vertices(mol)])
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
    push!(mol.gprops.descriptors.draw2d_bond_style, arr)
end


function coords2d(mol::SimpleMolGraph, i::Integer)
    dispatch_update!(mol)
    if !has_descriptor(mol, :coords2d)
        error("2D coordinates not available in this molecule type.")
    end
    cds = get_descriptor(mol, :coords2d)
    if length(cds) < i
        error("No available coords found. Run coordgen!(mol) to generate coords")
    elseif nv(mol) != length(cds[i])
        error("Coords size mismatch. Run coordgen!(mol) to regenerate coords")
    end
    return cds[i]
end
coords2d(mol::SimpleMolGraph) = coords2d(mol, 1)

function has_coords2d(mol::SimpleMolGraph)
    has_descriptor(mol, :coords2d) || return false
    isempty(get_descriptor(mol, :coords2d)) && return false
    length(get_descriptor(mol, :coords2d)[1]) == nv(mol) || return false
    return true
end


function coords3d(mol::SimpleMolGraph, i::Integer)
    dispatch_update!(mol)
    if !has_descriptor(mol, :coords3d)
        error("3D coordinates not available in this molecule type.")
    end
    cds = get_descriptor(mol, :coords3d)
    if length(cds) < i
        error("No available 3D coords found.")
    elseif nv(mol) != length(cds[i])
        error("Coords size mismatch.")
    end
    return cds[i]
end
coords3d(mol::SimpleMolGraph) = coords3d(mol, 1)

function has_coords3d(mol::SimpleMolGraph)
    has_descriptor(mol, :coords3d) || return false
    isempty(get_descriptor(mol, :coords3d)) && return false
    length(get_descriptor(mol, :coords3d)[1]) == nv(mol) || return false
    return true
end


function draw2d_bond_style(mol::SimpleMolGraph, i::Integer)
    dispatch_update!(mol)
    return get_descriptor(mol, :draw2d_bond_style)[i]
end
draw2d_bond_style(mol::SimpleMolGraph) = draw2d_bond_style(mol, 1)


function update_coords!(mol::SimpleMolGraph)
    # TODO: just remap nodes if still all existing vertices have coords.
    empty!(mol.gprops.descriptors.coords2d)
    empty!(mol.gprops.descriptors.coords3d)
    empty!(mol.gprops.descriptors.draw2d_bond_style)
    if hasfield(vproptype(mol), :coords)
        if !any(isnothing(get_prop(mol, i, :coords)) for i in vertices(mol))
            coords_from_sdf!(mol)
        end
    else
        coordgen!(mol)
    end
end



"""
    coordgen(mol::MolGraph) -> Tuple{Array{Int,1},Array{Int,1}}

Generate 2D coordinates by using Schrodinger's coordgenlibs.

This will returns a tuple of `coords` and `styles` arrays. `coords` is a size(n, 2) matrix
where n is atom count, which stores 2D coordinates (x, y) of each atoms.
`styles` is a size e vector of wedge notation of stereobond, where e is bond count.
"""
coordgen(mol::SimpleMolGraph) = coordgen(
    mol.graph, atom_number(mol), bond_order(mol),
    mol.gprops.stereocenter, mol.gprops.stereobond
)

function coordgen(mol::ReactiveMolGraph)
    dispatch_update!(mol)
    return coordgen(
        mol.graph, atom_number(mol), bond_order(mol),
        mol.gprops.stereocenter, mol.gprops.stereobond
    )
end

function coordgen!(mol::ReactiveMolGraph)
    # Initialize
    empty!(mol.gprops.descriptors.coords2d)
    empty!(mol.gprops.descriptors.draw2d_bond_style)
    # TODO: unspecified stereochem in SMILES
    coords, styles = coordgen(
        mol.graph, atom_number(mol), bond_order(mol),
        get_prop(mol, :stereocenter), get_prop(mol, :stereobond)
    )
    push!(mol.gprops.descriptors.coords2d, coords)
    push!(mol.gprops.descriptors.draw2d_bond_style, styles)
end

function coordgen(
        g::SimpleGraph{T}, atomnum::Vector{Int}, bondorder_::Vector{Int},
        stereocenters::Dict{T,Tuple{T,T,T,Bool}},
        stereobonds::Dict{Edge{T},Tuple{T,T,Bool}}) where T
    minmol = @ccall libcoordgen.getSketcherMinimizer()::Ptr{Cvoid}
    atoms = Ptr{Cvoid}[]
    bonds = Ptr{Cvoid}[]
    er = Dict(e => i for (i, e) in enumerate(edges(g)))

    # Atoms
    for a in atomnum
        a = a > 0 ? a : 9  # F as a placeholder of virtual atom
        atom = @ccall libcoordgen.setAtom(
            minmol::Ptr{Cvoid}, a::Cint)::Ptr{Cvoid}
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
