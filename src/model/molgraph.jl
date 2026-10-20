#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#




"""
    AbstractState

The base class of molecular model states.
"""
abstract type AbstractState end


"""
    MolState{T,F1,F2} <: AbstractState

The state container for `ReactiveMolGraph`.

`on_init::F1` and `on_update::F2` are callbacks for descriptor initialization and update.
Their types are part of the state type so that the callbacks are statically dispatched
(required for `juliac --trim`). The callbacks are fixed at construction; to replace them,
create a new molecule by `MolecularGraph.with_callbacks`.
"""
mutable struct MolState{T,F1,F2} <: AbstractState
    initialized::Bool
    has_updates::Bool
    has_new_edges::Bool  # as a SSSR recalculation flag
    disable_update::Bool  # ignore dispatch_update! (temporary set when running auto-preprocess)
    const on_init::F1
    const on_update::F2
end

function MolState{T}(
        ; initialized=false, has_updates=true, has_new_edges=true, disable_update=false,
        on_init::F1=default_on_init!, on_update::F2=default_on_update!
        ) where {T<:Integer,F1,F2}  # F1, F2: force specialization on callbacks
    return MolState{T,typeof(on_init),typeof(on_update)}(
        initialized, has_updates, has_new_edges, disable_update, on_init, on_update)
end

function Base.setproperty!(state::MolState, name::Symbol, x)
    name in (:on_init, :on_update) && error(
        "callbacks cannot be replaced in place. Use `MolecularGraph.with_callbacks(mol; $(name)=f)`")
    return setfield!(state, name, convert(fieldtype(typeof(state), name), x))
end

Base.copy(state::T) where T <: MolState = T(
    state.initialized, state.has_updates, state.has_new_edges, state.disable_update,
    state.on_init, state.on_update
)

# Default callbacks of a concrete state type (callbacks should be singleton functions)
function default_callbacks(::Type{MolState{T,F1,F2}}) where {T,F1,F2}
    (isdefined(F1, :instance) && isdefined(F2, :instance)) || error(
        "on_init and on_update should be given explicitly for $(MolState{T,F1,F2})")
    return (on_init=F1.instance, on_update=F2.instance)
end

StructUtils.structlike(::Type{<:MolState}) = false
StructUtils.lower(x::MolState) = Dict{String,Int}()  # state is not serialized
StructUtils.lift(::Type{<:MolState{T}}, x) where T = MolState{T}()


# Property update mechanisms

function get_descriptor(mol::ReactiveMolGraph, ::Val{F}) where F
    mol.state.disable_update || dispatch_update!(mol)
    return getfield(mol.gprops.descriptors, F)
end


"""
    remap!(gprop::SimpleMolProperty{T}, vmap::Vector{T}, edges::Vector{Edge{T}}) where T <: Integer -> Nothing

Remap vertices according to the given vmap (new_v = old_v -> vmap[old_v]).
"""
function remap!(gprop::SimpleMolProperty, args...)
    # vmap[old] -> new
    remap_fields!(gprop, args...)
    return
end

# Call `remap!(Val(field), container, args...)` for each field of the container.
# Unrolled at compile time so that the method for each field is statically dispatched.
@generated function remap_fields!(container, args...)
    calls = [:(remap!(Val($(QuoteNode(f))), container, args...)) for f in fieldnames(container)]
    return Expr(:block, calls..., :(return))
end

function remap!(::Val, gprop::SimpleMolProperty, args...)
    return
end

function remap!(::Val{:descriptors}, gprop::SimpleMolProperty, args...)
    remap_fields!(gprop.descriptors, args...)
    return
end



"""
    dispatch_update!(mol::ReactiveMolGraph) -> Nothing

Dispatch property auto-update mechanisms.

This is called by functions that returns calculated properties (descriptors).
If the vertices or edges has been added or removed, on_update callback is triggered
to recalculate descriptors, and then 2D wedges are updated to represent the stereocenters
([`assign_stereocenter_wedges!`](@ref)). If there are no changes, do nothing.
"""
function dispatch_update!(mol::ReactiveMolGraph)
    mol.state.has_updates || return
    mol.state.disable_update = true # to avoid infinite roop
    mol.state.on_update(mol)
    # keep 2D wedges consistent with the stereocenters (independent of the callback)
    assign_stereocenter_wedges!(mol)
    mol.state.disable_update = false
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

This mechanism is similar to `notify_updates!`, but is specific to the `add_edge!` and
`rem_edge!` methods.
Unlike other descriptors, the `sssr` descriptor does not require recalculation when vertices
are removed or added — remapping is sufficient in those cases.
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
    with_callbacks(mol::ReactiveMolGraph; on_init=nothing, on_update=nothing) -> ReactiveMolGraph

Return a molecule that shares the graph and properties with `mol` but has the given callbacks
(`nothing` keeps the current one). The state flags are copied and `on_init` is not called again.

Callbacks are part of the molecule type, so they cannot be replaced in place.
"""
function with_callbacks(mol::ReactiveMolGraph{T}; on_init::F1=nothing, on_update::F2=nothing
        ) where {T,F1,F2}  # F1, F2: force specialization on callbacks
    st = mol.state
    state = MolState{T}(;
        initialized=st.initialized, has_updates=st.has_updates,
        has_new_edges=st.has_new_edges, disable_update=st.disable_update,
        on_init=something(on_init, st.on_init), on_update=something(on_update, st.on_update))
    return _with_state(mol, state)
end


"""
    remap_gprops(mol::ReactiveMolGraph, vmap::Vector{T}, edges::Vector{Edge{T}}) -> MolProperty

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
    remap_gprops!(mol::ReactiveMolGraph, vmap::Vector{T}, edges::Vector{Edge{T}}) -> Nothing

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

function Base.setindex!(mol::ReactiveMolGraph, prop::AbstractElement, v::Integer)
    setindex!(mol.vprops, prop, v)
    notify_updates!(mol)
end

function Base.setindex!(mol::ReactiveMolGraph, prop::AbstractElement, e::Edge)
    setindex!(mol.eprops, prop, e)
    notify_updates!(mol)
end

function Base.setindex!(mol::ReactiveMolGraph{T,V,E}, prop::AbstractElement, u::T, v::T) where {T,V,E}
    setindex!(mol.eprops, prop, u_edge(T, u, v))
    notify_updates!(mol)
end


# old accessors (deprecated)
@deprecate set_prop!(mol::ReactiveMolGraph, v::Integer, value::AbstractElement) setindex!(mol, value, v)
@deprecate set_prop!(mol::ReactiveMolGraph, e::Edge, value::AbstractElement) setindex!(mol, value, e)


# 2D bond styles are indexed by the edge rank, so they are updated when an edge is added or
# removed (the new edge has no wedge)
function edge_styles_inserted!(mol::ReactiveMolGraph, e::Edge)
    has_descriptor(mol, Val(:draw2d_bond_style)) || return
    i = findfirst(==(e), collect(edges(mol.graph)))
    for styles in mol.gprops.descriptors.draw2d_bond_style
        length(styles) == ne(mol.graph) - 1 && insert!(styles.styles, i, :none)
    end
    return
end

function edge_styles_deleted!(mol::ReactiveMolGraph, i::Int)
    has_descriptor(mol, Val(:draw2d_bond_style)) || return
    for styles in mol.gprops.descriptors.draw2d_bond_style
        length(styles) == ne(mol.graph) + 1 && deleteat!(styles.styles, i)
    end
    return
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
    edge_styles_inserted!(mol, e)
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

If the edge is connected to a stereocenter or a stereobond (double bond), the stereochemistry
information is removed with a warning (see [`stereo_on_edge_removal!`](@ref)).
"""
function rem_u_edge!(mol::ReactiveMolGraph, e::Edge)
    # Can be directly called if src < dst is guaranteed.
    i = findfirst(==(e), collect(edges(mol.graph)))
    rem_edge!(mol.graph, e) || return false
    edge_styles_deleted!(mol, i)
    stereo_on_edge_removal!(mol, e)
    delete!(mol.eprops, e)
    notify_updates!(mol)
    notify_new_edges!(mol)  # rings may be broken
    return true
end


function Graphs.rem_vertex!(mol::ReactiveMolGraph, v::Integer)
    nv_ = nv(mol)
    has_vertex(mol.graph, v) || return false
    stereo_on_vertex_removal!(mol, [v])
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
    # remap is required even if the last vertex was removed (properties referring to it)
    vmap = collect(vertices(mol))
    nv_ == v || (vmap[v] = nv_)
    remap_gprops!(mol, vmap, old_edges)
    notify_updates!(mol)
    return true
end


function Graphs.rem_vertices!(mol::ReactiveMolGraph{T,V,E}, vs::Vector{T}) where {T,V,E}
    # TODO: if many vertices should be removed, induced_subgraph may be more efficient.
    stereo_on_vertex_removal!(mol, vs)
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
    # The subgraph has its own state (sharing the state with `mol` mixes up the update flags)
    state = copy(mol.state)
    state.has_updates = true  # descriptors of the subgraph should be recalculated
    # Rings may be broken in an edge-induced subgraph even if all of their vertices remain,
    # so SSSR should be recalculated.
    eltype(vlist_or_elist) <: Edge && (state.has_new_edges = true)
    return T(subg, vps, eps, newgp, state), vmap
end


function reactive_molgraph(
        g::SimpleGraph{T}, vprops::Dict{VertexKey{T},V},
        eprops::Dict{EdgeKey{T},E}, gprops::SimpleMolProperty{T}
        ; kwargs...) where {T,V,E}
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

reactive_molgraph(
    g::SimpleGraph{T}, vprops::Dict{T,V}, eprops::Dict{Edge{T},E}
    ; gprops=MolProperty{T}(), kwargs...
) where {T,V,E} = reactive_molgraph(
    g, convert(Dict{VertexKey{T},V}, vprops),
    convert(Dict{EdgeKey{T},E}, eprops), gprops; kwargs...
)

# from edge and property list (sdftomol, smilestomol interface)

function reactive_molgraph(
        edge_list::Vector{Edge{T}}, vprop_list::Vector{V}, eprop_list::Vector{E}
        ; gprops=MolProperty{T}(), kwargs...) where {T,V,E}
    g = SimpleGraph(edge_list)
    vps = Dict{VertexKey{T},V}(VertexKey(i) => v for (i, v) in enumerate(vprop_list))
    # eprop_list in edge_list order
    eps = Dict{EdgeKey{T},E}(EdgeKey(e) => eprop_list[i] for (i, e) in enumerate(edge_list))
    return reactive_molgraph(g, vps, eps, gprops; kwargs...)
end


"""
    MolGraph{T,V,E} <: ReactiveMolGraph{T,V,E}

Basic molecular graph type.

Usually `MolGraph`s are not directly constructed, but constructed by
molecule reader methods such as `smilestomol` or `sdfilereader`.
Atom and Bond properties of `MolGraph` are specific to the type of morecular reader
(e.g. sdfilereader generates `MolGraph{Int,SDFAtom,SDFBond}` instances).
"""
struct MolGraph{T<:Integer,V<:AbstractAtom,E<:AbstractBond,S<:MolState{T}} <: ReactiveMolGraph{T,V,E}
    graph::SimpleGraph{T}
    vprops::Dict{VertexKey{T},V}
    eprops::Dict{EdgeKey{T},E}
    gprops::MolProperty{T}
    state::S
end

# vprops and eprops are converted to Dict{VertexKey{T},V} and Dict{EdgeKey{T},E}
MolGraph{T,V,E}(
    g::SimpleGraph{T}, vprops::AbstractDict, eprops::AbstractDict,
    gprops::MolProperty{T}, state::S
) where {T,V,E,S<:MolState{T}} = MolGraph{T,V,E,S}(g, vprops, eprops, gprops, state)

function MolGraph{T,V,E}(args::Vararg{Any,N}; kwargs...) where {T,V,E,N}  # N: force specialization
    mol = MolGraph{T,V,E}(reactive_molgraph(args...; kwargs...)...)
    initialize!(mol)
    return mol
end

# Construct with the callbacks of the given concrete state type by default
MolGraph{T,V,E,S}(args::Vararg{Any,N}; kwargs...) where {T,V,E,S,N} = MolGraph{T,V,E}(
    args...; default_callbacks(S)..., kwargs...)::MolGraph{T,V,E,S}

_with_state(mol::MolGraph{T,V,E}, state::MolState{T}) where {T,V,E} = MolGraph{T,V,E}(
    mol.graph, mol.vprops, mol.eprops, mol.gprops, state)

MolGraph(
    g::SimpleGraph{T}, vprops::Dict{T,V}, eprops::Dict{Edge{T},E}; kwargs...
) where {T,V,E} = MolGraph{T,V,E}(g, vprops, eprops; kwargs...)
MolGraph(
    edge_list::Vector{Edge{T}}, vprop_list::Vector{V}, eprop_list::Vector{E}; kwargs...
) where {T,V,E} = MolGraph{T,V,E}(edge_list, vprop_list, eprop_list; kwargs...)

MolGraph{T,V,E}() where {T,V,E} = MolGraph(SimpleGraph{T}(), Dict{T,V}(), Dict{Edge{T},E}())
MolGraph{T,V,E,S}() where {T,V,E,S} = MolGraph{T,V,E,S}(SimpleGraph{T}(), Dict{T,V}(), Dict{Edge{T},E}())
MolGraph() = MolGraph{Int,SDFAtom,SDFBond}()


# MolGraph type aliases

"""
    SDFMolGraph

Alias of `MolGraph{Int,SDFAtom,SDFBond}`, the default molecule type of
[`sdftomol`](@ref) and [`sdfilereader`](@ref).
"""
const SDFMolGraph = MolGraph{Int,SDFAtom,SDFBond}
"""
    SMILESMolGraph

Alias of `MolGraph{Int,SMILESAtom,SMILESBond}`, the default molecule type of
[`smilestomol`](@ref).
"""
const SMILESMolGraph = MolGraph{Int,SMILESAtom,SMILESBond}
"""
    CommonChemMolGraph

Alias of `MolGraph{Int,CommonChemAtom,CommonChemBond}`, the molecule type for
the CommonChem JSON format.
"""
const CommonChemMolGraph = MolGraph{Int,CommonChemAtom,CommonChemBond}
