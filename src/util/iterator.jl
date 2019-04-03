#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    combinations, sortstablemax, sortstablemin


## combinations

struct Combinations{T}
    collection::Vector{T}
    indice::Vector{Int}
end

function Base.iterate(cmb::Combinations, state=nothing)
    lastidx = lastindex(cmb.collection)
    idcsize = length(cmb.indice)
    if cmb.indice[end] == lastidx
        if idcsize == 1
            return
        end
        n = 1
        while cmb.indice[end - n] == lastidx - n
            n += 1
            if n == idcsize
                return
            end
        end
        cmb.indice[idcsize - n] += 1
        for i in (idcsize - n + 1):idcsize
            cmb.indice[i] = cmb.indice[i - 1] + 1
        end
    else
        cmb.indice[end] += 1
    end
    return (cmb.collection[cmb.indice], state)
end

Base.eltype(::Type{Combinations{T}}) where {T} = Vector{T}
Base.IteratorSize(::Type{<:Combinations}) = Base.SizeUnknown()

"""
    combinations(iter, k::Int=2)

An iterator that yields `k`-combinations (card. k subsets of the collection)
"""
function combinations(iter, k::Int=2)
    col = collect(iter)
    T = eltype(col)
    if k < 0
        throw(DomainError(k, "k should not be negative"))
    elseif k > length(col)
        throw(ErrorException("$(k) is larger than the collection size"))
    elseif k == 0
        return (Int[],)
    else
        idc = collect(1:k)
        idc[end] -= 1
        return Combinations{T}(col, idc)
    end
end


function sortstablemax(iter; by=identity, kwargs...)
    isempty(iter) && return kwargs[:init]
    cmp(x, y) = ifelse(by(x) < by(y), y, x)
    return reduce(cmp, iter)
end


function sortstablemin(iter; by=identity, kwargs...)
    isempty(iter) && return kwargs[:init]
    cmp(x, y) = ifelse(by(x) > by(y), y, x)
    return reduce(cmp, iter)
end
