#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#



struct Stereocenter{T<:Integer}
    lookingFrom::T
    first::T
    second::T
    isclockwise::Bool
end

Base.Tuple(x::Stereocenter) = (x.lookingFrom, x.first, x.second, x.isclockwise)
Base.iterate(x::Stereocenter, state...) = iterate(Tuple(x), state...)
Base.:(==)(x::Stereocenter, y) = Tuple(x) == Tuple(y)
JSON.lower(x::Stereocenter) = collect(Tuple(x))


struct StereocenterMap{T<:Integer}
    mapping::Dict{VertexKey{T},Stereocenter{T}}
end

StereocenterMap{T}() where T = StereocenterMap{T}(Dict{VertexKey{T},Stereocenter{T}}())

function Base.iterate(x::StereocenterMap, state...)
    r = iterate(x.mapping, state...)
    return isnothing(r) ? nothing : (first(r[1]).key => last(r[1]), r[2])
end

Base.eltype(::StereocenterMap{T}) where T = Pair{T,Stereocenter{T}}
Base.:(==)(x::StereocenterMap, y::StereocenterMap) = x.mapping == y.mapping
Base.length(x::StereocenterMap) = length(x.mapping)
Base.copy(x::T) where T <: StereocenterMap = T(copy(x.mapping))
Base.keys(x::StereocenterMap) = getproperty.(keys(x.mapping), :key)
Base.haskey(x::StereocenterMap, k) = haskey(x.mapping, k)
Base.getindex(x::StereocenterMap{T}, k...) where T = getindex(x.mapping, VertexKey{T}.(k)...)
Base.setindex!(x::StereocenterMap{T}, v, k...) where T = setindex!(x.mapping, v, VertexKey{T}.(k)...)
Base.empty!(x::StereocenterMap) = empty!(x.mapping)
Base.merge!(x::T, y::T) where T <: StereocenterMap = merge!(x.mapping, y.mapping)

function remap(
        stereomap::StereocenterMap{T}, vmap::Vector{T}, edges::Vector{Edge{T}}) where T <: Integer
    revv = Dict(v => i for (i, v) in enumerate(vmap))
    newmap = StereocenterMap{T}()
    for (k, v) in stereomap
        isempty(setdiff([k, stereoneighbors(v)...], keys(revv))) || continue
        newmap[revv[k]] = Stereocenter{T}(revv[v.lookingFrom], revv[v.first], revv[v.second], v.isclockwise)
    end
    return newmap
end

function remap!(::Val{:stereocenter}, gprop::SimpleMolProperty{T}, args...) where T <: Integer
    newmap = remap(gprop.stereocenter, args...)
    empty!(gprop.stereocenter)
    merge!(gprop.stereocenter, newmap)
    return
end



struct Stereobond{T<:Integer}
    first::T
    second::T
    is_cis::Bool
end

Base.Tuple(x::Stereobond) = (x.first, x.second, x.is_cis)
Base.iterate(x::Stereobond, state...) = iterate(Tuple(x), state...)
Base.:(==)(x::Stereobond, y) = Tuple(x) == Tuple(y)
JSON.lower(x::Stereobond) = collect(Tuple(x))


struct StereobondMap{T<:Integer}
    mapping::Dict{EdgeKey{T},Stereobond{T}}
end

StereobondMap{T}() where T = StereobondMap{T}(Dict{EdgeKey{T},Stereobond{T}}())

function Base.iterate(x::StereobondMap, state...)
    r = iterate(x.mapping, state...)
    return isnothing(r) ? nothing : (first(r[1]).key => last(r[1]), r[2])
end

Base.:(==)(x::StereobondMap, y::StereobondMap) = x.mapping == y.mapping
Base.eltype(::StereobondMap{T}) where T = Pair{T,Stereobond{T}}
Base.length(x::StereobondMap) = length(x.mapping)
Base.copy(x::T) where T <: StereobondMap = T(copy(x.mapping))
Base.keys(x::StereobondMap) = getproperty.(keys(x.mapping), :key)
Base.haskey(x::StereobondMap, k) = haskey(x.mapping, k)
Base.getindex(x::StereobondMap{T}, k...) where T = getindex(x.mapping, EdgeKey{T}.(k)...)
Base.setindex!(x::StereobondMap{T}, v, k...) where T = setindex!(x.mapping, v, EdgeKey{T}.(k)...)
Base.empty!(x::StereobondMap) = empty!(x.mapping)
Base.merge!(x::T, y::T) where T <: StereobondMap = merge!(x.mapping, y.mapping)

function remap(
        stereomap::StereobondMap{T}, vmap::Vector{T}, edges::Vector{Edge{T}}) where T <: Integer
    revv = Dict(v => i for (i, v) in enumerate(vmap))
    newmap = StereobondMap{T}()
    for (k, v) in stereomap
        isempty(setdiff([src(k), dst(k), v.first, v.second], keys(revv))) || continue
        n1, n2 = revv[src(k)] < revv[dst(k)] ? (v.first, v.second) : (v.second, v.first)
        newmap[u_edge(T, revv[src(k)], revv[dst(k)])] = Stereobond(revv[n1], revv[n2], v.is_cis)
    end
    return newmap
end

function remap!(::Val{:stereobond}, gprop::SimpleMolProperty{T}, args...) where T <: Integer
    newmap = remap(gprop.stereobond, args...)
    empty!(gprop.stereobond)
    merge!(gprop.stereobond, newmap)
    return
end



struct Coords2d
    coords::Vector{Point2d}
end

Coords2d() = Coords2d([])
Coords2d(len::Int) = Coords2d(Vector(undef, len))
Base.:(==)(x::Coords2d, y::Coords2d) = x.coords == y.coords
Base.:(==)(x::Coords2d, y) = x.coords == y
Base.iterate(x::Coords2d, state...) = iterate(x.coords, state...)
Base.eltype(::Coords2d) = Point2d
Base.length(x::Coords2d) = length(x.coords)
Base.copy(x::Coords2d) = Coords2d(copy(x.coords))
Base.getindex(x::Coords2d, k...) = getindex(x.coords, k...)
Base.setindex!(x::Coords2d, v, k...) = setindex!(x.coords, v, k...)
Base.empty!(x::Coords2d) = empty!(x.coords)

StructUtils.structlike(::StructUtils.StructStyle, ::Type{Coords2d}) = false
JSON.lower(x::Coords2d) = [[p...] for p in x.coords]
JSON.lift(::Type{Coords2d}, x) = Coords2d([Point2d(p...) for p in x])

function remap(
        coords::Vector{Coords2d}, vmap::Vector{T}, edges::Vector{Edge{T}}) where T <: Integer
    revv = Dict(v => i for (i, v) in enumerate(vmap))
    container = Vector{Coords2d}(undef, length(coords))
    for (i, oldcds) in enumerate(coords)
        newcds = Coords2d(length(revv))
        for (oldp, newp) in revv
            newcds[newp] = oldcds[oldp]
        end
        container[i] = newcds
    end
    return container
end

function remap!(::Val{:coords2d}, gprop::SimpleMolProperty{T}, args...) where T <: Integer
    newmap = remap(gprop.coords2d, args...)
    empty!(gprop.coords2d)
    append!(gprop.coords2d, newmap)
    return
end



struct Coords3d
    coords::Vector{Point3d}
end

Coords3d() = Coords3d([])
Coords3d(len::Int) = Coords3d(Vector(undef, len))
Base.:(==)(x::Coords3d, y::Coords3d) = x.coords == y.coords
Base.:(==)(x::Coords3d, y) = x.coords == y
Base.iterate(x::Coords3d, state...) = iterate(x.coords, state...)
Base.eltype(::Coords3d) = Point3d
Base.length(x::Coords3d) = length(x.coords)
Base.copy(x::Coords3d) = Coords3d(copy(x.coords))
Base.getindex(x::Coords3d, k...) = getindex(x.coords, k...)
Base.setindex!(x::Coords3d, v, k...) = setindex!(x.coords, v, k...)
Base.empty!(x::Coords3d) = empty!(x.coords)

StructUtils.structlike(::StructUtils.StructStyle, ::Type{Coords3d}) = false
JSON.lower(x::Coords3d) = [[p...] for p in x.coords]
JSON.lift(::Type{Coords3d}, x) = Coords3d([Point3d(p...) for p in x])

function remap(
        coords::Vector{Coords3d}, vmap::Vector{T}, edges::Vector{Edge{T}}) where T <: Integer
    revv = Dict(v => i for (i, v) in enumerate(vmap))
    container = Vector{Coords3d}(undef, length(coords))
    for (i, oldcds) in enumerate(coords)
        newcds = Coords3d(length(revv))
        for (oldp, newp) in revv
            newcds[newp] = oldcds[oldp]
        end
        container[i] = newcds
    end
    return container
end

function remap!(::Val{:coords3d}, gprop::SimpleMolProperty{T}, args...) where T <: Integer
    newmap = remap(gprop.coords3d, args...)
    empty!(gprop.coords3d)
    append!(gprop.coords3d, newmap)
    return
end



struct Draw2dBondStyle
    styles::Vector{Symbol}
end

Draw2dBondStyle() = Draw2dBondStyle([])
Draw2dBondStyle(len::Int) = Draw2dBondStyle(Vector(undef, len))
Base.:(==)(x::Draw2dBondStyle, y::Draw2dBondStyle) = x.styles == y.styles
Base.:(==)(x::Draw2dBondStyle, y) = x.styles == y
Base.iterate(x::Draw2dBondStyle, state...) = iterate(x.styles, state...)
Base.eltype(::Draw2dBondStyle) = Symbol
Base.length(x::Draw2dBondStyle) = length(x.styles)
Base.copy(x::Draw2dBondStyle) = Draw2dBondStyle(copy(x.styles))
Base.getindex(x::Draw2dBondStyle, k...) = getindex(x.styles, k...)
Base.setindex!(x::Draw2dBondStyle, v, k...) = setindex!(x.styles, v, k...)
Base.empty!(x::Draw2dBondStyle) = empty!(x.styles)

StructUtils.structlike(::StructUtils.StructStyle, ::Type{Draw2dBondStyle}) = false
JSON.lower(x::Draw2dBondStyle) = string.(x.styles)
JSON.lift(::Type{Draw2dBondStyle}, x) = Draw2dBondStyle(Symbol.(x))

function remap(
        style::Vector{Draw2dBondStyle}, vmap::Vector{T}, edges::Vector{Edge{T}}) where T <: Integer
    container = Vector{Draw2dBondStyle}(undef, length(style))
    revv = Dict(v => i for (i, v) in enumerate(vmap))
    newedges = sort([u_edge(T, revv[src(e)], revv[dst(e)]) for e in edges
            if src(e) in vmap && dst(e) in vmap])
    olderank = Dict(e => i for (i, e) in enumerate(edges))
    emap = Dict(i => edge_rank(olderank, vmap[src(e)], vmap[dst(e)]) for (i, e) in enumerate(newedges))
    inv = Dict{Symbol,Symbol}(:up => :revup, :revup => :up, :down => :revdown, :revdown => :down)
    for (i, styles) in enumerate(style)
        cont = Draw2dBondStyle(length(emap))
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
    return container
end

function remap!(::Val{:draw2d_bond_style}, gprop::SimpleMolProperty, args...)
    newmap = remap(gprop.draw2d_bond_style, args...)
    empty!(gprop.draw2d_bond_style)
    append!(gprop.draw2d_bond_style, newmap)
    return
end



# Common interfaces of AbstractProperty

# Metadata specific shorthands (e.g. mol["compound_id"] = "CP000001")

Base.getindex(mol::SimpleMolGraph, key::String) = mol.gprops.metadata[key]

function Base.setindex!(mol::SimpleMolGraph, value::String, key::String)
    # skip dispatch (Metadata update would not affect graph state)
    mol.gprops.metadata[key] = value
end

# old accessors (deprecated)
get_prop(mol::SimpleMolGraph, key::String) = mol[key]
has_prop(mol::SimpleMolGraph, key::String) = haskey(mol.gprops.metadata, key)
set_prop!(mol::SimpleMolGraph, key::String, value::String) = setindex!(mol, value, key)


# Descriptors

"""
    MolDescriptor{T} <: SimpleMolProperty{T}

Container of calculated secondary properties (descriptor) compatible with `ReactiveMolGraph`.
"""
mutable struct MolDescriptor{T} <: SimpleMolProperty{T}
    # standardized atom charges and bond orders
    atom_charge::Vector{Int}
    bond_order::Vector{Int}
    # cached descriptors that are expensive or frequently used
    sssr::Vector{Vector{T}}
    lone_pair::Vector{Int}
    apparent_valence::Vector{Int}
    valence::Vector{Int}
    is_ring_aromatic::Vector{Bool}
    # coordinates
    coords2d::Vector{Coords2d}
    coords3d::Vector{Coords3d}
    draw2d_bond_style::Vector{Draw2dBondStyle}  # wedge notation in drawing
end

function MolDescriptor{T}(;
        atom_charge::Vector{Int}=Int[],
        bond_order::Vector{Int}=Int[],
        sssr::Vector{Vector{T}}=Vector{T}[],
        lone_pair::Vector{Int}=Int[],
        apparent_valence::Vector{Int}=Int[],
        valence::Vector{Int}=Int[],
        is_ring_aromatic::Vector{Bool}=Bool[],
        coords2d::Vector{Coords2d}=Coords2d[],
        coords3d::Vector{Coords3d}=Coords3d[],
        draw2d_bond_style::Vector{Draw2dBondStyle}=Draw2dBondStyle[]) where T <: Integer
    return MolDescriptor{T}(
        atom_charge, bond_order, sssr, lone_pair, apparent_valence,
        valence, is_ring_aromatic, coords2d, coords3d, draw2d_bond_style
    )
end

Base.copy(desc::T) where T <: MolDescriptor = T(
    copy(desc.atom_charge), copy(desc.bond_order), copy_vec_of_vec(desc.sssr),
    copy(desc.lone_pair), copy(desc.apparent_valence),
    copy(desc.valence), copy(desc.is_ring_aromatic),
    [copy(c) for c in desc.coords2d], [copy(c) for c in desc.coords3d],
    [copy(c) for c in desc.draw2d_bond_style]
)



# Properties

"""
    MolProperty{T} <: SimpleMolProperty{T}

Container of graph-level molecule properties compatible with `ReactiveMolGraph`.
"""
struct MolProperty{T} <: SimpleMolProperty{T}
    stereocenter::StereocenterMap{T}
    stereobond::StereobondMap{T}
    pyrrole_like::Vector{T}  # pyrrole H position for SMILES kekulization
    smarts_input::String  # TODO: should be metadata
    smarts_lexical_succ::Vector{Vector{T}}  # lexical index used for stereochem
    descriptors::MolDescriptor{T}
    metadata::OrderedDict{String,String}  # e.g. SDFile option fields
    logs::Dict{String,String}  # Parse errors
end

function MolProperty{T}(;
        stereocenter::StereocenterMap{T}=StereocenterMap{T}(),
        stereobond::StereobondMap{T}=StereobondMap{T}(),
        pyrrole_like::Vector{T}=T[],
        smarts_input::String="",
        smarts_lexical_succ::Vector{Vector{T}}=Vector{T}[],
        descriptors::MolDescriptor{T}=MolDescriptor{T}(),
        metadata::OrderedDict{String,String}=OrderedDict{String,String}(),
        logs::Dict{String,String}=Dict{String,String}()) where T <: Integer
    return MolProperty{T}(
        stereocenter, stereobond, pyrrole_like, smarts_input, smarts_lexical_succ,
        descriptors, metadata, logs
    )
end


Base.copy(prop::T) where T <: MolProperty = T(
    copy(prop.stereocenter), copy(prop.stereobond), copy(prop.pyrrole_like),
    prop.smarts_input, copy_vec_of_vec(prop.smarts_lexical_succ), copy(prop.descriptors),
    copy(prop.metadata), copy(prop.logs)
)
