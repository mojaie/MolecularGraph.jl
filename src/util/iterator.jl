#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    combinations


import Base: iterate, eltype, IteratorSize


struct Combinations{T}
    collection::Vector{T}
    indice::Vector{Int}
end

function Combinations(iter, num::Int)
    col = collect(iter)
    T = eltype(col)
    if num < 0 || num > length(col)
        throw(ValueError("invalid value $(num)"))
    elseif num == 0
        return (Int[],)
    else
        idc = collect(1:num)
        idc[end] -= 1
        return Combinations{T}(col, idc)
    end
end

"""
    combinations(iter, k::Int)

An iterator that yields `k`-combinations (subsets of the `iter` elements)
"""
combinations(iter, k::Int=2) = Combinations(iter, k)

function iterate(cmb::Combinations, state=nothing)
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

eltype(::Type{Combinations{T}}) where {T} = Vector{T}
IteratorSize(::Type{<:Combinations}) = Base.SizeUnknown()
