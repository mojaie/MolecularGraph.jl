#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

"""
    ReactiveMolGraph{T<:Integer,V<:AbstractElement,E<:AbstractElement} <: SimpleMolGraph{T}

The base class of molecule model which have auto-update mechanism of properties.

Typically `ReactiveMolGraph` should have the following properties:
- `graph`: molecular graph topology in `Graphs.SimpleGraph`.
- `vprops`: `Vector` of atom properties (e.g. SDFAtom, SMILESAtom)
- `eprops`: `Vector` of bond properties (e.g. SDFBond, SMILESBond)
- `gprops`: graph-level properties and stored descriptors (e.g. stereocenter)
- `states`: update flags and callback functions for `reactive` property update

"""
abstract type ReactiveMolGraph{T<:Integer,V<:AbstractElement,E<:AbstractElement} <: SimpleMolGraph{T} end

Base.:(==)(g::ReactiveMolGraph, h::ReactiveMolGraph
    ) = g.graph == h.graph && g.vprops == h.vprops && g.eprops == h.eprops && g.gprops == h.gprops

vproptype(::Type{<:ReactiveMolGraph{T,V,E}}) where {T,V,E} = V
eproptype(::Type{<:ReactiveMolGraph{T,V,E}}) where {T,V,E} = E

Base.copy(mol::T) where T <: ReactiveMolGraph = T(
    copy(mol.graph), copy(mol.vprops), copy(mol.eprops), copy(mol.gprops), copy(mol.state))


"""
    MolState{T,F1,F2} <: AbstractState

The state container for `ReactiveMolGraph`.
"""
mutable struct MolState{T,F1,F2} <: AbstractState
    initialized::Bool
    has_updates::Bool
    has_new_edges::Bool  # as a SSSR recalculation flag
    force_calculate::Bool  # ignore cached descriptors (temporary set when running auto-preprocess)
    on_init::F1
    on_update::F2
end


function MolState{T}(
        ; initialized=false, has_updates=true, has_new_edges=true, force_calculate=false,
        on_init=default_on_init!, on_update=default_on_update!) where T <: Integer
    return MolState{T,typeof(on_init),typeof(on_update)}(
        initialized, has_updates, has_new_edges, force_calculate, on_init, on_update)
end


Base.copy(state::T) where T <: MolState = T(
    state.initialized, state.has_updates, state.has_new_edges, state.force_calculate,
    state.on_init, state.on_update
)


# Property update mechanisms

function get_descriptor(mol::ReactiveMolGraph, field::Symbol)
    dispatch_update!(mol)
    return getproperty(mol[:descriptors], field)
end

function has_descriptor(mol::ReactiveMolGraph, field::Symbol)
    mol.state.force_calculate && return false
    return hasproperty(mol[:descriptors], field)
end


"""
    remap!(container::SimpleMolProperty{T}, vmap::DVector{T}, edges::Vector{Edge{T}}) where T <: Integer

Remap vertices according to the given vmap (new_v = old_v -> vmap[old_v]).
"""
function remap!(container::SimpleMolProperty{T}, vmap::Vector{T},
        edges::Vector{Edge{T}}) where T <: Integer
    # vmap[old] -> new
    for sym in fieldnames(typeof(container))
        remap!(Val(sym), container, vmap, edges)
    end
    return
end

function remap!(::Val, container::SimpleMolProperty{T}, vmap::Vector{T},
        edges::Vector{Edge{T}}) where T <: Integer
    return
end

function remap!(::Val{:descriptors}, gprop::SimpleMolProperty{T},
        vmap::Vector{T}, edges::Vector{Edge{T}}) where T <: Integer
    for sym in fieldnames(typeof(gprop.descriptors))
        remap!(Val(sym), gprop.descriptors, vmap, edges)
    end
    return
end


"""
    dispatch_update!(mol::ReactiveMolGraph) -> Nothing

Dispatch property auto-update mechanisms.

This is called by functions that returns calculated properties (descriptors).
If the vertices or edges has been added or removed, on_update callback is triggered
to recalculate descriptors. If there are no changes, just return the stored descriptors.
"""
function dispatch_update!(mol::ReactiveMolGraph)
    mol.state.has_updates || return
    mol.state.force_calculate = true # to avoid infinite roop
    mol.state.on_update(mol)
    mol.state.force_calculate = false
    mol.state.has_updates = false  
    mol.state.has_new_edges = false
    return
end


"""
    notify_updates!(mol::ReactiveMolGraph) -> Nothing

Set :has_updates flag to inform property updaters that recalculation is required.

This is called by graph manipulation functions (add, remove or set).
"""
function notify_updates!(mol::ReactiveMolGraph)
    # TODO: flag to inactivate auto-update
    mol.state.has_updates = true
    return
end


"""
    notify_new_edges!(mol::ReactiveMolGraph) -> Nothing

Set `:has_new_edges` flag to inform property updaters that recalculation is required.

This mechanism is similar to `notify_updates!`, but is specific to the `add_edge!` method.  
Unlike other descriptors, the `sssr` descriptor does not require recalculation when vertices or edge
are removed, or when new vertices are added — remapping is sufficient in those cases.
Therefore, the `:has_new_edges` flag is used to avoid unnecessary and costly minimum cycle recalculation.
"""
function notify_new_edges!(mol::ReactiveMolGraph)
    # TODO: flag to inactivate auto-update
    mol.state.has_new_edges = true
    return
end


"""
    initialize!(mol::ReactiveMolGraph) -> Nothing

Initialize `ReactiveMolGraph`

Custom `on_init` function is called inside this method.
"""
function initialize!(mol::ReactiveMolGraph)
    mol.state.initialized || mol.state.on_init(mol)
    mol.state.initialized = true
    dispatch_update!(mol)  # check on-update callback errors
    return
end


"""
    remap_gprops(mol::ReactiveMolGraph, vmap::Dict{T,T}) -> MolProperty

Return updated `MolProperty`.

Only `Graph.rem_vertex!` variants should call this function to remap vertices
after removal.
"""
function remap_gprops(mol::ReactiveMolGraph, vmap::Vector{T},
        edges::Vector{Edge{T}}) where T <: Integer
    gprop = copy(mol.gprops)
    remap!(gprop, vmap, edges)
    return gprop
end

"""
    remap_gprops!(mol::ReactiveMolGraph, vmap::Dict{T,T}) -> Nothing

Update `MolProperty` in place.

Only `Graph.rem_vertex!` variants should call this function to remap vertices
after removal.
"""
function remap_gprops!(mol::ReactiveMolGraph, vmap::Vector{T},
        edges::Vector{Edge{T}}) where T <: Integer
    remap!(mol.gprops, vmap, edges)
    return
end


"""
    default_on_init!(mol::ReactiveMolGraph) -> Nothing

Default `on_init` callback for `ReactiveMolGraph`

This should be implemented indivisually to each molecule readers.
"""
function default_on_init!(mol::ReactiveMolGraph)
    return  # No initialization by default
end


"""
    default_on_update!(mol::ReactiveMolGraph) -> Nothing

Default `on_update` callback for `ReactiveMolGraph`

This should be implemented indivisually to each molecule readers.
"""
function default_on_update!(mol::ReactiveMolGraph)
    return  # No auto-update callbacks by default
end


# Property accessors

"""
    set_prop!(mol::ReactiveMolGraph, v::Integer, value::AbstractAtom) -> Nothing
    set_prop!(mol::ReactiveMolGraph, e::Edge, value::AbstractBond) -> Nothing

Set properties (vertex or edge attributes).
"""
function set_prop!(mol::ReactiveMolGraph, v::Integer, value::AbstractElement)
    mol.vprops[v] = value
    notify_updates!(mol)
end

function set_prop!(mol::ReactiveMolGraph, e::Edge, value::AbstractElement)
    mol.eprops[e] = value
    notify_updates!(mol)
end


"""
    add_u_edge!(mol::ReactiveMolGraph, e::Edge, prop::AbstractElement) -> Bool

Add undirected edge to the `ReactiveMolGraph` and notifies the update to updator functions.

This is called internally by the API function `add_edge!(mol, e, prop)`.
`add_edge!(mol, e, prop)` always add the sorted edge to the `undirected` molecular graph.
`add_u_edge!` can be directly called if src < dst is guaranteed.
"""
function add_u_edge!(mol::ReactiveMolGraph, e::Edge, prop::AbstractElement)
    # Can be directly called if src < dst is guaranteed.
    add_edge!(mol.graph, e) || return false
    mol.eprops[e] = prop
    notify_updates!(mol)
    notify_new_edges!(mol)
    return true
end


function Graphs.add_vertex!(mol::ReactiveMolGraph{T,V,E}, prop::V) where {T,V,E}
    add_vertex!(mol.graph) || return false
    mol.vprops[nv(mol.graph)] = prop
    notify_updates!(mol)
    return true
end


"""
    rem_u_edge!(mol::ReactiveMolGraph, e::Edge) -> Bool

Remove undirected edge to the `ReactiveMolGraph` and notifies the update to updator functions.

This is called internally by the API function `rem_edge!(mol, e)`.
`rem_edge!(mol, e)` always remove the sorted edge from the `undirected` molecular graph.
`rem_u_edge!` can be directly called if src < dst is guaranteed.
"""
function rem_u_edge!(mol::ReactiveMolGraph, e::Edge)
    # Can be directly called if src < dst is guaranteed.
    rem_edge!(mol.graph, e) || return false
    delete!(mol.eprops, e)
    notify_updates!(mol)
    return true
end


function Graphs.rem_vertex!(mol::ReactiveMolGraph, v::Integer)
    nv_ = nv(mol)
    old_edges = collect(edges(mol))
    rem_vertex!(mol.graph, v) || return false
    # last index node is re-indexed to the removed node
    nv_ == v || begin mol.vprops[v] = mol.vprops[nv_] end
    delete!(mol.vprops, nv_)
    for e in old_edges
        (src(e) == v || dst(e) == v) && delete!(mol.eprops, e)
	    ((src(e) == v && dst(e) == nv_) || (src(e) == nv_ && dst(e) == v)) && continue
        nv_ == v && continue
        src(e) == nv_ && begin mol.eprops[u_edge(mol, v, dst(e))] = mol.eprops[e]; delete!(mol.eprops, e) end
        dst(e) == nv_ && begin mol.eprops[u_edge(mol, src(e), v)] = mol.eprops[e]; delete!(mol.eprops, e) end
    end
    if nv_ != v
        vmap = collect(vertices(mol))
        vmap[v] = nv_
        remap_gprops!(mol, vmap, old_edges)
    end
    notify_updates!(mol)
    return true
end


function Graphs.rem_vertices!(mol::ReactiveMolGraph{T,V,E}, vs::Vector{T}) where {T,V,E}
    # TODO: if many vertices should be removed, induced_subgraph may be more efficient.
    old_edges = collect(edges(mol))
    vmap = rem_vertices!(mol.graph, vs)
    # remap vertex properties
    for (i, v) in enumerate(vmap)
        i == v && continue
        mol.vprops[i] = mol.vprops[v]
    end
    for v in (nv(mol)+1):(nv(mol)+length(vs))
        delete!(mol.vprops, v)
    end
    # remap edge properties
    new_edges = collect(edges(mol))
    for e in new_edges
        _e = u_edge(T, vmap[src(e)], vmap[dst(e)])
        e == _e && continue
        mol.eprops[e] = mol.eprops[_e]
    end
    for e in keys(mol.eprops)
        e in new_edges || delete!(mol.eprops, e)
    end
    remap_gprops!(mol, vmap, old_edges)
    notify_updates!(mol)
    return vmap
end


function _induced_subgraph(mol::T, vlist_or_elist) where {T<:ReactiveMolGraph}
    # In many cases, induced_subgraph(mol.graph, vlist_or_elist) is sufficient
    subg, vmap = induced_subgraph(mol.graph, vlist_or_elist)
    vps = Dict{eltype(T),vproptype(T)}()
    for v in vertices(subg)
        vps[v] = mol.vprops[vmap[v]]
    end
    eps = Dict{Edge{eltype(T)},eproptype(T)}()
    for e in edges(subg)
        eps[e] = mol.eprops[u_edge(mol, vmap[src(e)], vmap[dst(e)])]
    end
    newgp = remap_gprops(mol, vmap, collect(edges(mol)))
    notify_updates!(mol)
    return T(subg, vps, eps, newgp, mol.state), vmap
end


function reactive_molgraph(
        g::SimpleGraph{T}, vprops::Dict{T,V}, eprops::Dict{Edge{T},E};
        gprops=MolProperty{T}(), kwargs...) where {T,V,E}
    if nv(g) > length(vprops)
        error("Mismatch in the number of nodes and node properties")
    elseif ne(g) != length(eprops)
        error("Mismatch in the number of edges and edge properties")
    end
    # expand fadjlist for vprops of isolated nodes
    for _ in nv(g):(length(vprops) - 1)
        push!(g.fadjlist, T[])
    end
    config = MolState{T}(;kwargs...)
    return (g, vprops, eprops, gprops, config)
end

# from edge and property list (sdftomol, smilestomol interface)

function reactive_molgraph(
        edge_list::Vector{Edge{T}}, vprop_list::Vector{V}, eprop_list::Vector{E}
        ; kwargs...) where {T,V,E}
    g = SimpleGraph(edge_list)
    vps = Dict{T,V}(i => v for (i, v) in enumerate(vprop_list))
    # eprop_list in edge_list order
    eps = Dict{Edge{T},E}(e => eprop_list[i] for (i, e) in enumerate(edge_list))
    return reactive_molgraph(g, vps, eps; kwargs...)
end


"""
    MolGraph{T,V,E} <: ReactiveMolGraph{T,V,E}

Basic molecular graph type.

Usually `MolGraph`s are not directly constructed, but constructed by
molecule reader methods such as `smilestomol` or `sdfilereader`.
Atom and Bond properties of `MolGraph` are specific to the type of morecular reader
(e.g. sdfilereader generates `MolGraph{Int,SDFAtom,SDFBond}` instances).
"""
struct MolGraph{T<:Integer,V<:AbstractAtom,E<:AbstractBond} <: ReactiveMolGraph{T,V,E}
    graph::SimpleGraph{T}
    vprops::Dict{T,V}
    eprops::Dict{Edge{T},E}
    gprops::MolProperty{T}
    state::MolState{T}
end

function MolGraph{T,V,E}(args...; kwargs...) where {T,V,E}
    mol = MolGraph{T,V,E}(reactive_molgraph(args...; kwargs...)...)
    initialize!(mol)
    return mol
end

MolGraph(
    g::SimpleGraph{T}, vprops::Dict{T,V}, eprops::Dict{Edge{T},E}; kwargs...
) where {T,V,E} = MolGraph{T,V,E}(g, vprops, eprops; kwargs...)
MolGraph(
    edge_list::Vector{Edge{T}}, vprop_list::Vector{V}, eprop_list::Vector{E}; kwargs...
) where {T,V,E} = MolGraph{T,V,E}(edge_list, vprop_list, eprop_list; kwargs...)

MolGraph{T,V,E}() where {T,V,E} = MolGraph(SimpleGraph{T}(), Dict{T,V}(), Dict{Edge{T},E}())
MolGraph() = MolGraph{Int,SDFAtom,SDFBond}()


# MolGraph type aliases

const SDFMolGraph = MolGraph{Int,SDFAtom,SDFBond}
const SMILESMolGraph = MolGraph{Int,SMILESAtom,SMILESBond}
const CommonChemMolGraph = MolGraph{Int,CommonChemAtom,CommonChemBond}
