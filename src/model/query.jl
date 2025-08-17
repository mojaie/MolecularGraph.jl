#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#


"""
    QueryNode

Query components

- query literals (arg -> key[arg] == value).
- logical operators (arg -> q1[arg] && q2[arg]).
- tautology query (arg -> true/false).
"""
struct QueryNode <: AbstractQueryNode
    operator::Symbol  # :and, :or, :not, :any, :eq, (:gt, :lt)
    key::Symbol  # keys for :eq, (:gt, :lt)
    value::String  # :any(Bool), :eq, (:gt, :lt)
end

function QueryNode(;
        operator::Union{AbstractString,Symbol}=:eq, 
        key::Union{AbstractString,Symbol}=:none,
        value::AbstractString="")
    return QueryNode(Symbol(operator), Symbol(key), value)
end

QueryNode(d::Dict{String,Any}
    ) = QueryNode(; NamedTuple((Symbol(k), v) for (k, v) in d)...)
QueryNode(d::Dict{Symbol,Any}
    ) = QueryNode(; NamedTuple((k, v) for (k, v) in d)...)

qand() = QueryNode(:and, :none, "")
qor() = QueryNode(:or, :none, "")
qnot() = QueryNode(:not, :none, "")
qanytrue() = QueryNode(:any, :none, "")
qeq(key, value) = QueryNode(:eq, key, value)
qtrue(key) = qeq(key, "true")


function Base.isless(q::QueryNode, r::QueryNode)
    # Node sorting required by truth table generator
    q.operator < r.operator && return true
    q.operator > r.operator && return false
    q.key < r.key && return true
    q.key > r.key && return false
    return q.value < r.value
end

function to_dict(::Val{:default}, q::QueryNode)
    rcd = Dict{String,Any}()
    q.operator === :eq || setindex!(rcd, string(q.operator), "operator")
    q.key === :none || setindex!(rcd, string(q.key), "key")
    q.value == "" || setindex!(rcd, q.value, "value")
    return rcd
end


function querytree(edges::Vector{Tuple{T,T}}, props::Vector{U}) where {T,U}
    g = SimpleDiGraph(Edge{T}.(edges))
    vprops = Dict{T,U}(i => p for (i, p) in enumerate(props))
    # expand fadjlist and badjlist for vprops of isolated nodes
    for _ in nv(g):(length(vprops) - 1)
        push!(g.fadjlist, T[])
        push!(g.badjlist, T[])
    end
    return g, vprops
end


function querytree(::Type{T}, ::Type{U}, data::Dict{String,Any}) where {T,U}
    edges = NTuple{2,T}[(e...,) for e in data["edges"]]
    props = U.(data["vprops"])
    return querytree(edges, props)
end


function to_dict(fmt::Val{:default}, qtree::QueryTree)
    return Dict{String,Any}(
        "edges" => [[src(e), dst(e)] for e in edges(qtree.graph)],
        "vprops" => [to_dict(fmt, qtree.vprops[i]) for i in vertices(qtree.graph)]
    )
end


struct QueryAtom{T<:Integer,U<:QueryNode} <: QueryTree{T,U}
    graph::SimpleDiGraph{T}
    vprops::Dict{T,U}
end

struct QueryBond{T<:Integer,U<:QueryNode} <: QueryTree{T,U}
    graph::SimpleDiGraph{T}
    vprops::Dict{T,U}
end

Base.copy(qtree::T) where T <: QueryAtom = T(copy(qtree.graph), copy(qtree.vprops))
Base.copy(qtree::T) where T <: QueryBond = T(copy(qtree.graph), copy(qtree.vprops))

QueryAtom{T,U}(edges::Vector, props::Vector
    ) where {T<:Integer,U<:QueryNode} = QueryAtom{T,U}(querytree(edges, props)...)
QueryBond{T,U}(edges::Vector, props::Vector
    ) where {T<:Integer,U<:QueryNode} = QueryBond{T,U}(querytree(edges, props)...)
QueryAtom(edges::Vector, props::Vector) = QueryAtom(querytree(edges, props)...)
QueryBond(edges::Vector, props::Vector) = QueryBond(querytree(edges, props)...)
QueryAtom() = QueryAtom(Tuple{Int,Int}[], QueryNode[])
QueryBond() = QueryBond(Tuple{Int,Int}[], QueryNode[])
QueryAtom{T,U}(data::Dict{String,Any}
    ) where {T<:Integer,U<:QueryNode} = QueryAtom{T,U}(querytree(T, U, data)...)
QueryBond{T,U}(data::Dict{String,Any}
    ) where {T<:Integer,U<:QueryNode} = QueryBond{T,U}(querytree(T, U, data)...)
QueryAtom(data::Dict{String,Any}) = QueryAtom(querytree(Int, QueryNode, data)...)
QueryBond(data::Dict{String,Any}) = QueryBond(querytree(Int, QueryNode, data)...)



@kwdef mutable struct QueryMolDescriptor{T} <: SimpleMolProperty{T}
    coords2d::Vector{Vector{Point2d}} = Vector{Point2d}[]
    coords3d::Vector{Vector{Point3d}} = Vector{Point3d}[]
end

Base.copy(desc::T) where T <: QueryMolDescriptor = T(
    copy_vec_of_vec(desc.coords2d), copy_vec_of_vec(desc.coords3d)
)


@kwdef struct QueryMolProperty{T} <: SimpleMolProperty{T}
    stereocenter::Dict{T,Tuple{T,T,T,Bool}} = Dict{T,Tuple{T,T,T,Bool}}()
    stereobond::Dict{Edge{T},Tuple{T,T,Bool}} = Dict{Edge{T},Tuple{T,T,Bool}}()
    smarts_input::String = ""
    smarts_lexical_succ::Vector{Vector{T}} = Vector{T}[]  # lexical index used for stereochem
    smarts_connectivity::Vector{Vector{T}} = Vector{T}[]  # SMARTS connectivity query
    descriptors::QueryMolDescriptor{T} = QueryMolDescriptor{T}()
    # Graph-level metadata properties (e.g. SDFile option fields)
    metadata::OrderedDict{String,String} = OrderedDict{String,String}()
    # Parse errors
    logs::Dict{String,String} = Dict{String,String}()
end


Base.copy(prop::T) where T <: QueryMolProperty = T(
    copy(prop.stereocenter), copy(prop.stereobond), prop.smarts_input,
    copy_vec_of_vec(prop.smarts_lexical_succ), copy_vec_of_vec(prop.smarts_connectivity),
    copy(prop.descriptors), copy(prop.metadata), copy(prop.logs)
)


reconstruct(::Val{:descriptors}, ::Type{QueryMolProperty{T}}, @nospecialize(data)
    ) where T = reconstruct(QueryMolDescriptor{T}, data)


"""
    QueryMolGraph{T,V,E} <: ReactiveMolGraph{T,V,E}

Basic molecular graph type.
"""
struct QueryMolGraph{T<:Integer,V<:QueryTree,E<:QueryTree} <: ReactiveMolGraph{T,V,E}
    graph::SimpleGraph{T}
    vprops::Dict{T,V}
    eprops::Dict{Edge{T},E}
    gprops::QueryMolProperty{T}
    state::MolState{T}
end

function QueryMolGraph{T,V,E}(args...; kwargs...) where {T,V,E}
    mol = QueryMolGraph{T,V,E}(
        reactive_molgraph(args...; gprops=QueryMolProperty{T}(), kwargs...)...)
    initialize!(mol)
    return mol
end

QueryMolGraph(
    g::SimpleGraph{T}, vprops::Dict{T,V}, eprops::Dict{Edge{T},E}; kwargs...
) where {T,V,E} = QueryMolGraph{T,V,E}(g, vprops, eprops; kwargs...)
QueryMolGraph(
    edge_list::Vector{Edge{T}}, vprop_list::Vector{V}, eprop_list::Vector{E}; kwargs...
) where {T,V,E} = QueryMolGraph{T,V,E}(edge_list, vprop_list, eprop_list; kwargs...)

QueryMolGraph{T,V,E}() where {T,V,E} = QueryMolGraph(SimpleGraph{T}(), Dict{T,V}(), Dict{Edge{T},E}())
QueryMolGraph() = QueryMolGraph{Int,QueryAtom,QueryBond}()


"""
    canonical(qtree::QueryTree) -> Tuple{String,Vector{Int}}

Return cannonical representation of String and the vertex order.
"""
function canonical(qtree::QueryTree)
    root(qtree)  # connectivity check
    cnn = Vector{String}(undef, nv(qtree.graph))
    order = Vector{Int}[[i] for i in 1:nv(qtree.graph)]
    revtopo = topological_sort(reverse(qtree.graph))
    for v in revtopo
        succs = outneighbors(qtree.graph, v)
        vp = qtree.vprops[v]
        if isempty(succs)
            cnn[v] = string(vp.operator, vp.key, vp.value, "/")
        else
            o = sortperm(cnn[succs])
            concat = reduce(*, cnn[succs[o]])
            cnn[v] = string(vp.operator, vp.key, vp.value, concat)
            append!(order[v], order[succs[o]]...)
        end
    end
    return cnn[revtopo[end]], order[revtopo[end]]
end


function Base.:(==)(q::QueryTree, r::QueryTree)
    typeof(q) !== typeof(r) && return false
    nv(q.graph) != nv(r.graph) && return false
    nv(q.graph) == 0 && return true
    ne(q.graph) != ne(r.graph) && return false
    return canonical(q)[1] == canonical(r)[1]
end


function Base.hash(qtree::QueryTree, h::UInt)
    nv(qtree.graph) == 0 && return hash([])
    hashes = Vector{UInt}(undef, nv(qtree.graph))
    revtopo = topological_sort(reverse(qtree.graph))
    for v in revtopo
        succs = outneighbors(qtree.graph, v)
        if isempty(succs)
            hashes[v] = hash(qtree.vprops[v])
        else
            hs = reduce(hash, sort(hashes[succs]))
            hashes[v] = hash(qtree.vprops[v], hs)
        end
    end
    return hash(typeof(qtree), hashes[revtopo[end]])
end


"""
    root(qtree::QueryTree) -> Integer

Return root node of the tree.
"""
root(qtree::QueryTree) = only(findall(iszero, indegree(qtree.graph)))


function add_qnode!(qtree::QueryTree, prop::QueryNode)
    add_vertex!(qtree.graph)
    qtree.vprops[nv(qtree.graph)] = prop
    return nv(qtree.graph)
end

function add_qnode!(qtree::QueryTree, src::Integer, prop::QueryNode)
    add_vertex!(qtree.graph)
    qtree.vprops[nv(qtree.graph)] = prop
    newnode = nv(qtree.graph)
    add_edge!(qtree.graph, src, newnode)
    return newnode
end


function set_qnode!(qtree::QueryTree, i::Integer, prop::QueryNode)
    qtree.vprops[i] = prop
    return
end


function rem_qnode!(qtree::QueryTree, v::Integer)
    nv_ = nv(qtree.graph)
    rem_vertex!(qtree.graph, v)
    # last index node is re-indexed to the removed node
    nv_ == v || setindex!(qtree.vprops, qtree.vprops[nv_], v)
    delete!(qtree.vprops, nv_)
    return
end


function rem_qnodes!(qtree::QueryTree, vs::Vector{<:Integer})
    vmap = rem_vertices!(qtree.graph, vs)
    for (i, v) in enumerate(vmap)
        i == v && continue
        qtree.vprops[i] = qtree.vprops[v]
    end
    for v in (nv(qtree.graph)+1):(nv(qtree.graph)+length(vs))
        delete!(qtree.vprops, v)
    end
    return vmap
end


add_qedge!(qtree::QueryTree, src::Integer, dst::Integer) = add_edge!(qtree.graph, src, dst)
rem_qedge!(qtree::QueryTree, src::Integer, dst::Integer) = rem_edge!(qtree.graph, src, dst)


"""
    querypropmap(qtree::QueryTree) -> Dict{Symbol,Vector{QueryNode}}

Parse QueryLiteral tree and put QueryLiterals into bins labeled with their literal keys.
"""
function querypropmap(qtree::QueryTree)
    pmap = Dict{Symbol,Vector{QueryNode}}()
    for (i, p) in qtree.vprops
        p.operator === :eq || continue
        if !haskey(pmap, p.key)
            pmap[p.key] = QueryNode[]
        end
        push!(pmap[p.key], p)
    end
    return pmap
end


"""
    generate_queryfunc(tree, props) -> Function

Generate query truthtable function from QueryLiteral tree and the property vector.

The query truthtable function take a Vector{Bool} of length equal to `props` and
returns output in Bool.
"""
function generate_queryfunc(qtree::QueryTree, props::Vector{QueryNode})
    nv(qtree.graph) == 0 && return arr -> false
    func = Dict{Int,Function}()
    revtopo = topological_sort(reverse(qtree.graph))
    for v in revtopo
        p = qtree.vprops[v]
        succs = outneighbors(qtree.graph, v)
        if isempty(succs)
            if p.operator === :eq
                idx = findfirst(x -> x == p, props)
                func[v] = arr -> arr[idx]
            elseif p.operator === :any
                func[v] = arr -> true
            else
                @assert false p.operator
            end
            continue
        end
        sfunc = [func[s] for s in succs]
        if p.operator === :not
            func[v] = arr -> ~only(sfunc)(arr)
        elseif p.operator === :and
            func[v] = arr -> all(f(arr) for f in sfunc)
        elseif p.operator === :or
            func[v] = arr -> any(f(arr) for f in sfunc)
        else
            @assert false
        end
    end
    return func[revtopo[end]]
end


"""
    specialize_nonaromatic!(q::QueryMolGraph) -> Nothing

Convert `[#atomnumber]` queries connected to explicit single bonds to be non-aromatic
(e.g. -[#6]- -> -C-).

Should be applied before `remove_hydrogens!`.
This function is intended for generalization of PAINS query in PubChem dataset.
"""
function specialize_nonaromatic!(qmol::QueryMolGraph{T,V,E}) where {T,V,E}
    aromsyms = Set([:B, :C, :N, :O, :P, :S, :As, :Se])
    exqs = Set(E(
            [(1, 2), (1, 3), (3, 4)],
            [qand(), qeq(:order, string(i)), qnot(), qtrue(:isaromatic)]
        ) for i in 1:3)
    exbonds = Dict{Edge{Int},Int}()
    for e in edges(qmol)
        qmol[e] in exqs || continue
        exbonds[e] = parse(Int, qmol[e].vprops[2].value)
    end
    for i in vertices(qmol)
        qtree = qmol[i]
        qonly = qtree.vprops[1]
        # e.g. [#6], [#7], [#8]...
        (nv(qtree.graph) == 1 && qonly.key === :symbol) || continue
        # number of explicitly non-aromatic incident bonds
        cnt = sum(get(exbonds, u_edge(qmol, i, nbr), 0) for nbr in neighbors(qmol, i); init=0)
        if qonly.value === "C"
            cnt -= 1  # carbon allows one non-aromatic
        end
        if Symbol(qonly.value) in aromsyms && cnt > 0
            set_prop!(qmol, i, V(
                [(1, 2), (1, 3), (3, 4)],
                [qand(), qonly, qnot(), qtrue(:isaromatic)]
            ))
        end
    end
    return
end


"""
    resolve_not_hydrogen! -> Nothing

Return the molecular query with hydrogen nodes removed.

This function is intended for generalization of PAINS query in PubChem dataset.
"""
function resolve_not_hydrogen!(qtree::QueryTree{T,U}) where {T,U}
    toremove = T[]
    for (i, prop) in qtree.vprops
        prop.operator === :not || continue
        succ = only(outneighbors(qtree.graph, i))
        sprop = qtree.vprops[succ]
        if sprop.key === :symbol && sprop.value === "H"
            push!(toremove, succ)
            set_qnode!(qtree, i, qanytrue())
        end
    end
    rem_qnodes!(qtree, toremove)
    return
end


"""
    remove_hydrogens!(q::MolGraph) -> Nothing

Remove hydrogens from the molecular query.

Should be applied after `specialize_nonaromatic!`.
This function is intended for generalization of PAINS query in PubChem dataset.
"""
function remove_hydrogens!(qmol::QueryMolGraph{T,V,E}) where {T,V,E}
    # count H nodes and mark H nodes to remove
    hnodes = T[]
    hcntarr = zeros(Int, nv(qmol))
    for n in vertices(qmol)
        qtree = qmol[n]
        resolve_not_hydrogen!(qtree)  # [!#1] -> [*]
        nv(qtree.graph) == 1 || continue
        (qtree.vprops[1].key === :symbol && qtree.vprops[1].value === "H") || continue
        hcntarr[only(neighbors(qmol, n))] += 1
        push!(hnodes, n)
    end
    for n in setdiff(vertices(qmol), hnodes)  # heavy atom nodes
        # e.g. C([H])([H]) is the same as [C;!H0;!H1]
        hs = collect(0:(hcntarr[n] - 1))
        isempty(hs) && continue
        qtree = qmol[n]
        rt = root(qtree)
        node = add_qnode!(qtree, qand())
        add_qedge!(qtree, node, rt)
        for i in hs
            c = add_qnode!(qtree, node, qnot())
            add_qnode!(qtree, c, qeq(:total_hydrogens, string(i)))
        end
        set_prop!(qmol, n, qtree)
    end
    rem_vertices!(qmol, hnodes)
    return
end


"""
    preprocess!(qmol::QueryMolGraph, smarts::String)

The default SMARTS preprocessor called after when SMARTS is parsed.
"""
function preprocess!(qmol::QueryMolGraph)
    smarts = qmol.gprops.smarts_input
    if occursin(r"-", smarts) && occursin(r"#[1-9]", smarts)
        specialize_nonaromatic!(qmol)
    end
    if occursin(r"#1", smarts)
        remove_hydrogens!(qmol)
    end
    return
end



# Not used

"""
    optimize_query!(tree::QueryTree{T,V}) where {T<:Integer,V<:QueryNode}) -> Nothing

Optimize query tree.

- Absorption
  - (:and => (A, A)) -> A
  - (:or => (A, A, B)) -> (:or => (A, B))
  - (:and => (:any => true, A)) -> A
  - :or => (:any => true, A) -> :any => true
- `:not` has the highest precedence in SMARTS, but only in the case like [!C],
  De Morgan's law will be applied to remove `:and` under `:not`.
  (e.g. :not => (:and => (:atomsymbol => :C, :isaromatic => false)
   -> :or => (:not => (:atomsymbol => :C), isaromatic => true)
"""
function optimize!(tree::QueryTree{T,U}) where {T,U}
    # TODO: possibly unnecessary

    # Remove `:and` under `:not`
    toremove = []
    # Traverse from leafs to the root
    revtopo = topological_sort(reverse(tree.graph))
    for i in revtopo
        tree.vprops[i].operator === :not || continue
        c = only(outneighbors(tree.graph, i))
        tree.vprops[c].operator === :and || continue
        succs = outneighbors(tree.graph, c)
        set_qnode!(tree, i, qor())
        push!(toremove, c)
        for succ in succs  # by-pass qor -> qnot -> successors
            newc = add_qnode!(tree, i, qnot())
            add_qedge!(tree, newc, succ)
        end
    end
    # Absorption
    for i in revtopo
        prop = tree.vprops[i]
        succs = outneighbors(tree.graph, i)
        if prop.operator in (:eq, :any)
            @assert isempty(succs)
            continue
        elseif prop.operator === :not
            @assert length(succs) == 1
            continue
        end
        @assert prop.operator in (:and, :or) && !isempty(succs)
        smap = Dict{T,U}()
        anys = T[]
        for s in succs
            tree.vprops[s] in values(smap) && continue
            if tree.vprops[s].operator === :any
                push!(anys, s)
            else
                smap[s] = tree.vprops[s]
            end
        end
        parent = inneighbors(tree.graph, i)
        append!(toremove, setdiff(succs, keys(smap)))  # remove duplicates
        if length(smap) == 1
            # :and => (:foo,:any) -> :foo
            # :and/or => (:foo,) -> :foo
            # :and/or => (:foo,:foo) -> :foo
            push!(toremove, i)
            isempty(parent) || add_qedge!(tree, parent, only(keys(smap)))
        elseif (isempty(smap) || prop.operator === :or) && !isempty(anys)
            # :and => (:any,:any) -> :any
            # :or => (:foo,:any) -> :any
            push!(toremove, i)
            isempty(parent) || add_qedge!(tree, parent, qanytrue())
        end
    end
    rem_qnodes!(tree, toremove)
    return
end
