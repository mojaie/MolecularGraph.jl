#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

const MIN_TIMEOUT = 10  # minimum value of clique detection timeout
get_expire(timeout::Real) = UInt64(time_ns() + timeout * 1_000_000_000)


# Cliques

abstract type AbstractCliquesState end

abstract type CliquesState{T<:Integer} <: AbstractCliquesState end
abstract type ConnCliquesState{T<:Integer} <: AbstractCliquesState end


function max_clique_postproc!(state::AbstractCliquesState)
    if length(state.Q) > length(state.maxsofar)
        state.maxsofar = copy(state.Q)
    end
end

function cliques_postproc!(state::AbstractCliquesState)
    push!(state.cliques, copy(state.Q))
end


mutable struct AllCliquesState{T,U,F} <: CliquesState{T}
    graph::SimpleGraph{T}
    targetsize::Union{Int,Nothing} # if Q reached the size, finish and return cliques found so far.
    expire::Union{UInt64,Nothing} # UInt64, nanoseconds
    postprocess::F
    Q::Vector{T}
    cliques::Vector{U}  # Vector{T} (MCIS) or Dict{Edge{T},Edge{T}} (MCES)
    status::Symbol
end

AllCliquesState{T,U}(g::SimpleGraph{T};
    timeout=MIN_TIMEOUT, targetsize=nothing, postprocess=cliques_postproc!, kwargs...
) where {T,U} = AllCliquesState{T,U,typeof(postprocess)}(
    g, targetsize, get_expire(timeout), postprocess, [], U[], :ongoing)


mutable struct MaxCliqueState{T,U,F} <: CliquesState{T}
    graph::SimpleGraph{T}
    targetsize::Union{Int,Nothing} # if Q reached the size, finish and return cliques found so far.
    expire::Union{UInt64,Nothing} # UInt64, nanoseconds
    postprocess::F
    Q::Vector{T}
    maxsofar::U
    status::Symbol
end

MaxCliqueState{T,U}(g::SimpleGraph{T};
    timeout=MIN_TIMEOUT, targetsize=nothing, postprocess=max_clique_postproc!, kwargs...
) where {T,U} = MaxCliqueState{T,U,typeof(postprocess)}(
    g, targetsize, get_expire(timeout), postprocess, [], U(), :ongoing)


function expand!(state::CliquesState{T}, subg::Set{T}, cand::Set{T}) where T
    (state.status == :timedout || state.status == :targetreached) && return
    if isempty(subg)
        # Report max clique
        state.postprocess(state)
        return
    elseif state.expire !== nothing && time_ns() > state.expire
        state.status = :timedout
        return
    elseif state.targetsize !== nothing && length(state.Q) >= state.targetsize
        state.status = :targetreached
        state.postprocess(state)
        return
    end
    # pivot = getpivot(state.graph, subg, cand)
    # maybe better to avoid expensive pivot selection in typical MCS problem.
    pivot = sortstablemax(subg, by=x->degree(state.graph, x))
    copv = setdiff(cand, neighbors(state.graph, pivot))
    for q in copv
        push!(state.Q, q)
        qnbrs = neighbors(state.graph, q)
        candq = intersect(cand, qnbrs)
        if (state isa AllCliquesState
                || length(state.Q) + length(candq) > length(state.maxsofar))
            subgq = intersect(subg, qnbrs)
            expand!(state, subgq, candq)
        end
        pop!(cand, q)
        pop!(state.Q)
    end
end

function getpivot(g::CliquesState{T}, subg::Set{T}, cand::Set{T}) where T
    maxadj = -1
    pv = 0
    for s in sort(collect(subg), by=x->degree(g, x), rev=true)
        degree(g, s) <= maxadj && break
        adjcnt = length(intersect(cand, neighbors(g, s)))
        if adjcnt > maxadj
            maxadj = adjcnt
            pv = s
        end
    end
    return pv
end


"""
    all_maximal_cliques(::Type{U}, g::SimpleGraph{T}; kwargs...
        ) where {T,U} -> (Vector{U}, Symbol)

Calculate maximal cliques.

# Reference

1. Tomita, E., Tanaka, A., & Takahashi, H. (2006). The worst-case time
   complexity for generating all maximal cliques and computational experiments.
   Theoretical Computer Science, 363(1), 28–42.
   https://doi.org/10.1016/J.TCS.2006.06.015
1. Cazals, F., & Karande, C. (2008). A note on the problem of reporting maximal
   cliques. Theoretical Computer Science, 407(1–3), 564–568.
   https://doi.org/10.1016/j.tcs.2008.05.010
"""
function all_maximal_cliques(::Type{U}, g::SimpleGraph{T}; kwargs...) where {T,U}
    state = AllCliquesState{T,U}(g; kwargs...)
    expand!(state, Set(vertices(g)), Set(vertices(g)))
    if state.status == :ongoing
        state.status = :done
    end
    return state.cliques, state.status
end
all_maximal_cliques(g::SimpleGraph{T}; kwargs...
    ) where T = all_maximal_cliques(Vector{T}, g; kwargs...)


"""
    maximum_clique(::Type{U}, g::SimpleGraph{T}; kwargs...) where {T,U} -> (U, Symbol)

Calculate maximum clique.
"""
function maximum_clique(::Type{U}, g::SimpleGraph{T}; kwargs...) where {T,U}
    state = MaxCliqueState{T,U}(g; kwargs...)
    expand!(state, Set(vertices(g)), Set(vertices(g)))
    if state.status == :ongoing
        state.status = :done
    end
    return state.maxsofar, state.status
end
maximum_clique(g::SimpleGraph{T}; kwargs...
    ) where T = maximum_clique(Vector{T}, g; kwargs...)



# Connected cliques (c-cliques)

function connfuncgen(modprod::SimpleGraph{T}, eattr::Dict{Edge{T},Bool}) where T <: Integer
    # connectivity adjlist
    conn = [Set{T}() for _ in vertices(modprod)]
    disconn = [Set{T}() for _ in vertices(modprod)]
    for i in vertices(modprod)
        for nbr in neighbors(modprod, i)
            container = eattr[u_edge(T, i, nbr)] ? conn : disconn
            push!(container[i], nbr)
        end
    end
    return n -> (conn[n], disconn[n])
end



mutable struct AllConnCliquesState{T,U,F1,F2} <: ConnCliquesState{T}
    graph::SimpleGraph{T}
    targetsize::Union{Int,Nothing} # if Q reached the size, finish and return cliques found so far.
    expire::Union{UInt64,Nothing} # UInt64, nanoseconds
    connfunc::F1
    postprocess::F2
    Q::Vector{T}
    cliques::Vector{U}
    status::Symbol
end

AllConnCliquesState{T,U,F}(g::SimpleGraph{T}, connfunc::F;
    timeout=MIN_TIMEOUT, targetsize=nothing, postprocess=cliques_postproc!, kwargs...
) where {T,U,F} = AllConnCliquesState{T,U,F,typeof(postprocess)}(g, targetsize,
    get_expire(timeout), connfunc, postprocess, [], U[], :ongoing)


mutable struct MaxConnCliqueState{T,U,F1,F2} <: ConnCliquesState{T}
    graph::SimpleGraph{T}
    targetsize::Union{Int,Nothing} # if Q reached the size, finish and return cliques found so far.
    expire::Union{UInt64,Nothing} # UInt64, nanoseconds
    connfunc::F1
    postprocess::F2
    Q::Vector{T}
    maxsofar::U
    status::Symbol
end

MaxConnCliqueState{T,U,F}(g::SimpleGraph{T}, connfunc::F;
    timeout=MIN_TIMEOUT, targetsize=nothing, postprocess=max_clique_postproc!, kwargs...
) where {T,U,F} = MaxConnCliqueState{T,U,F,typeof(postprocess)}(g, targetsize,
    get_expire(timeout), connfunc, postprocess, T[], U(), :ongoing)



function expand!(
        state::ConnCliquesState{T},
        P::Set{T}, Q::Set{T}, X::Set{T}, Y::Set{T}) where T
    (state.status == :timedout || state.status == :targetreached) && return
    if isempty(P) && isempty(X)
        # Report max clique
        state.postprocess(state)
        return
    elseif state.expire !== nothing && time_ns() > state.expire
        state.status = :timedout
        return
    elseif state.targetsize !== nothing && length(state.Q) >= state.targetsize
        state.status = :targetreached
        state.postprocess(state)
        return
    end
    pivot = sortstablemax(union(P, X), by=x->degree(state.graph, x))
    copv = setdiff(P, neighbors(state.graph, pivot))
    while !isempty(copv)
        n = pop!(copv)
        push!(state.Q, n)  # TODO: state.Q is R in orignal literature
        nbrs = neighbors(state.graph, n)
        conn::Set{T}, disconn::Set{T} = state.connfunc(n)
        Qnew = intersect(Q, disconn)
        Pnew = union(intersect(P, nbrs), intersect(Q, conn))
        Ynew = intersect(Y, disconn)
        Xnew = union(intersect(X, nbrs), intersect(Y, conn))
        if (state isa AllConnCliquesState
                || length(state.Q) + length(Pnew) + length(Qnew) > length(state.maxsofar))
            expand!(state, Pnew, Qnew, Xnew, Ynew)
        end
        push!(X, n)
        pop!(state.Q)
    end
end


function init_conn_cliques!(state::ConnCliquesState{T}) where T
    nodes = Set(vertices(state.graph))
    done = T[]
    for n in vertices(state.graph)
        push!(state.Q, n)  # TODO: Q_ should be renamed
        conn::Set{T}, disconn::Set{T} = state.connfunc(n)
        P = intersect(setdiff(nodes, done), conn)
        Q_ = intersect(setdiff(nodes, done), disconn)
        X = intersect(conn, done)
        Y = intersect(disconn, done)
        expand!(state, P, Q_, X, Y)
        push!(done, n)
        pop!(state.Q)
    end
    if state.status == :ongoing
        state.status = :done
    end
end


"""
    all_maximal_conn_cliques(::Type{U}, g::SimpleGraph{T}, eattrs::Dict;
        kwargs...) where {T,U} -> (Vector{U}, Symbol)

Calculate maximal connected cliques.

# Reference

1. Cazals, F., & Karande, C. (2005). An algorithm for reporting maximal
   c-cliques. Theoretical Computer Science, 349(3), 484–490.
   https://doi.org/10.1016/j.tcs.2005.09.038

"""
function all_maximal_conn_cliques(::Type{U}, g::SimpleGraph{T};
        connfunc=n->(Set{T}(), Set(neighbors(g, n))), kwargs...) where {T,U}
    state = AllConnCliquesState{T,U,typeof(connfunc)}(g, connfunc; kwargs...)
    init_conn_cliques!(state)
    return state.cliques, state.status
end
all_maximal_conn_cliques(g::SimpleGraph{T}; kwargs...
    ) where T = all_maximal_conn_cliques(Vector{T}, g; kwargs...)


"""
    maximum_conn_clique(::Type{U}, g::SimpleGraph{T}, eattrs::Dict;
        kwargs...) where {T,U} -> (U, Symbol)

Calculate maximal connected cliques.
"""
function maximum_conn_clique(::Type{U}, g::SimpleGraph{T};
    connfunc=n->(Set{T}(), Set(neighbors(g, n))), kwargs...) where {T,U}
    state = MaxConnCliqueState{T,U,typeof(connfunc)}(g, connfunc; kwargs...)
    init_conn_cliques!(state)
    return state.maxsofar, state.status
end
maximum_conn_clique(g::SimpleGraph{T}; kwargs...
    ) where T = maximum_conn_clique(Vector{T}, g; kwargs...)




# Approximate cliques (Experimental)


function ramsey_R2(g::SimpleGraph{T}, subg, pivot) where T
    # TODO: refactor
    isempty(subg) && return T[], T[]
    vset = setdiff(Set(subg), pivot)
    nbrs = collect(intersect(vset, neighbors(g, pivot)))
    nnbrs = collect(setdiff(vset, nbrs))
    c_1, i_1 = isempty(nbrs) ? (T[], T[]) : ramsey_R2(g, nbrs, nbrs[1])
    c_2, i_2 = isempty(nnbrs) ? (T[], T[]) : ramsey_R2(g, nnbrs, nnbrs[1])
    push!(c_1, pivot)
    push!(i_2, pivot)
    return argmax(length, [c_1, c_2]), argmax(length, [i_1, i_2])
end


"""
    approx_maximum_clique(g) -> Vector{Int}

Return approximate maximal cliques.

Reference:
- Boppana, R., & Halldórsson, M. M. (1992).
- https://networkx.org/documentation/stable/_modules/networkx/algorithms/approximation/clique.html#max_clique

"""
function approx_maximum_clique(g::SimpleGraph{T}, root=1) where T
    # TODO: not so performant, exact algorithm with targetsize or timeout option looks better
    nv(g) == 0 && return T[]
    g_ = complement(g)
    c_i, i_i = ramsey_R2(g_, collect(vertices(g_)), root)
    isets = [i_i]
    vmap_all = collect(vertices(g_))
    while nv(g_) > 0
        vmap = rem_vertices!(g_, c_i)
        vmap_all = vmap_all[vmap]
        c_i, i_i = ramsey_R2(g_, collect(vertices(g_)), 1)
        isempty(i_i) || push!(isets, vmap_all[i_i])
    end
    return sortstablemax(isets, by=length, init=[])
end
