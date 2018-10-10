#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    Ring,
    topology!,
    minifyring!

import Base: ==, getindex, size


mutable struct Ring <: AbstractVector{UInt16}
    arr::AbstractVector{UInt16}
    function Ring(nodes::AbstractVector)
        # Canonicalize
        (fst, fstidx) = findmin(nodes)
        succidx = fstidx == lastindex(nodes) ? 1 : fstidx + 1
        succ = nodes[succidx]
        predidx = fstidx == 1 ? lastindex(nodes) : fstidx - 1
        pred = nodes[predidx]
        cp = succ < pred ? copy(nodes) : reverse(copy(nodes))
        while cp[1] != fst
            pushfirst!(cp, pop!(cp))
        end
        new(cp)
    end
    Ring() = new()
end

size(r::Ring) = length(r.arr)
getindex(r::Ring, i::Int) = r.arr[i]
==(r::Ring, s::Ring) = r.arr == s.arr


function topology!(mol::MolecularGraph)
    nodeset = Set(a.index for a in atomvector(mol))
    biconnected = Dict()
    isolated = []
    while !isempty(nodeset)
        start = pop!(nodeset)
        stack = [start]
        pred = Dict(start => start)
        used = Dict(start => Set())
        root = Dict(start => start)
        while !isempty(stack)
            tail = pop!(stack)
            for nbr in keys(neighbors(mol, tail))
                if nbr ∉ keys(used) # New node
                    pred[nbr] = tail
                    push!(stack, nbr)
                    used[nbr] = Set(tail)
                    root[nbr] = nbr
                elseif nbr in stack # Cycle found
                    pn = used[nbr]
                    cyc = [nbr, tail]
                    p = pred[tail]
                    ed = pred[nbr]
                    root[nbr] = root[tail] = root[ed]
                    while p ∉ pn # Backtrack
                        push!(cyc, p)
                        root[p] = root[ed]
                        if p in keys(biconnected) # Append scaffold to new cycle
                            if root[ed] ∉ keys(biconnected)
                                biconnected[root[ed]] = []
                            end
                            append!(biconnected[root[ed]], biconnected[p])
                            delete!(biconnected, p)
                        end
                        p = pred[p]
                    end
                    push!(cyc, p)
                    if root[ed] ∉ keys(biconnected) # Append new cycle to scaffold
                        biconnected[root[ed]] = []
                    end
                    push!(biconnected[root[ed]], cyc)
                    push!(used[nbr], tail)
                end
            end
        end
        push!(isolated, collect(keys(pred)))
        # print(pred)
        setdiff!(nodeset, keys(pred))
    end
    mol.rings = Vector{Ring}()
    mol.scaffolds = Vector{Vector{UInt16}}()
    for cycles in values(biconnected)
        rcnt = length(mol.rings)
        append!(mol.rings, [Ring(c) for c in cycles])
        push!(mol.scaffolds, collect(rcnt + 1 : rcnt + length(cycles)))
    end
    mol.isolated = Vector{Vector{UInt16}}(sort(isolated, by=(x -> -length(x)))[2:end])
    push!(mol.descriptors, "Topology")
    return
end


function minifyring!(mol::MolecularGraph; verbose=false)
    required_descriptor(mol, "Topology")
    for cyc_idx in mol.scaffolds
        rings = sort([mol.rings[c] for c in cyc_idx], by=length)
        minified = []
        cnt = 0
        while !isempty(rings)
            cnt += 1
            if cnt > 100
                push!(mol.descriptors, "MinifiedRing")
                throw(OperationError("Ring minimization failed"))
            end
            r = popfirst!(rings)
            init_r = r
            if verbose
                print("$(length(r)) Ring:$(r)\n")
            end
            for m in minified
                if verbose
                    print("$(length(m)) Minified:$(m)\n")
                end
                resolved = resolve_inclusion(r, m)
                if resolved != nothing
                    if verbose
                        print(
                            "$(length(resolved[1])), $(length(resolved[2]))",
                            " Resolved:$(resolved[1]) $(resolved[2])\n"
                        )
                    end
                    r = resolved[1]
                end
            end
            if verbose
                print("$(length(r)) New ring:$(r)\n")
            end
            if length(r) == length(init_r) # no further minimization
                push!(minified, r)
            else
                push!(rings, r)
            end
        end
        for c in cyc_idx
            mol.rings[c] = pop!(minified)
        end
    end
    push!(mol.descriptors, "MinifiedRing")
end


function resolve_inclusion(a::Ring, b::Ring)
    rev = length(a.arr) > length(b.arr)
    (lt, bg) = rev ? (copy(b.arr), copy(a.arr)) : (copy(a.arr), copy(b.arr))
    isec = intersect(lt, bg)
    isize = length(isec)
    if isize == 3 && length(lt) == 4
        # TODO: cubane special case
    elseif isize != length(lt) && isize <= length(lt) / 2 + 1
        return nothing # Already minimum
    end
    while bg[1] ∉ isec || bg[end] in isec
        pushfirst!(bg, pop!(bg))
    end
    while lt[1] != bg[1]
        pushfirst!(lt, pop!(lt))
    end
    if bg[2] != lt[2] # Reverse
        lt = cat(lt[1], reverse(lt[2:end]), dims=1)
    end
    if bg[1:isize] != lt[1:isize]
        return nothing # Aboid irregular minification
    end
    bx = bg[isize:end]
    lx = lt[isize + 1 : end]
    ca = cat(bx, bg[1], reverse(lx), dims=1)
    new_bg = Ring(ca)
    rev ? (new_bg, Ring(lt)) : (Ring(lt), new_bg)
end
