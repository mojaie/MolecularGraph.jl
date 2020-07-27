#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    combinations, sortstablemax, sortstablemin


## combinations

struct Combinations
    n::Int
    k::Int
    indice::Vector{Int}
end

function Base.iterate(cmb::Combinations, state=nothing)
    if cmb.indice[end] == cmb.n
        cmb.k == 1 && return
        i = 1
        while cmb.indice[end - i] == cmb.n - i
            i += 1
            i == cmb.k && return
        end
        cmb.indice[cmb.k - i] += 1
        for j in (cmb.k - i + 1):cmb.k
            cmb.indice[j] = cmb.indice[j - 1] + 1
        end
    else
        cmb.indice[end] += 1
    end
    return (cmb.indice, state)
end

Base.eltype(::Type{<:Combinations}) = Vector{Int}
Base.IteratorSize(::Type{<:Combinations}) = Base.SizeUnknown()


"""
    combinations(n::Int, k::Int=2)

An iterator that yields ``k``-combinations (cardinality ``k`` subsets of the cardinality ``n`` set) as an array of integers.
"""
function combinations(n::Int, k::Int=2)
    if k < 0
        throw(DomainError(k, "k should not be negative"))
    elseif k > n
        throw(ErrorException(
            "$(k) should be equal or less than the collection size"))
    elseif k == 0
        return (Int[],)
    else
        idc = collect(1:k)
        idc[end] -= 1
        return Combinations(n, k, idc)
    end
end



function sortstablemax(iter; by=identity, kwargs...)
    cmp(x, y) = by(x) < by(y) ? y : x
    return reduce(cmp, iter; kwargs...)
end


function sortstablemin(iter; by=identity, kwargs...)
    cmp(x, y) = by(x) > by(y) ? y : x
    return reduce(cmp, iter; kwargs...)
end
