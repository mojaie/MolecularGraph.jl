#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    MolGraphTopology,
    molgraph_topology,
    molgraph_topology!


mutable struct MolGraphTopology <: Descriptor
    cycles::Vector
    bicomps::Vector
    connectivity::Vector
    cyclemap::Dict
end


function molgraph_topology(mol::MolecularGraph)
    nodeset = Set(a.index for a in atomvector(mol))
    candidates = Dict()
    connectivity = []
    while !isempty(nodeset)
        start = pop!(nodeset)
        stack = [start]
        pred = Dict(start => start)
        used = Dict(start => Set())
        root = Dict(start => start)
        while !isempty(stack)
            tail = pop!(stack)
            for nbr in keys(neighbors(mol, tail))
                if !(nbr in keys(used))
                    # New node
                    pred[nbr] = tail
                    push!(stack, nbr)
                    used[nbr] = Set(tail)
                    root[nbr] = nbr
                elseif nbr in stack
                    # Cycle found
                    pn = used[nbr]
                    cyc = [nbr, tail]
                    p = pred[tail]
                    ed = pred[nbr]
                    root[nbr] = root[tail] = root[ed]
                    while !(p in pn)
                        # Backtrack
                        push!(cyc, p)
                        root[p] = root[ed]
                        if p in keys(candidates)
                            # Append bicomp to new cycle
                            if !(root[ed] in keys(candidates))
                                candidates[root[ed]] = []
                            end
                            append!(candidates[root[ed]], candidates[p])
                            delete!(candidates, p)
                        end
                        p = pred[p]
                    end
                    push!(cyc, p)
                    if !(root[ed] in keys(candidates))
                        # Append new cycle to bicomp
                        candidates[root[ed]] = []
                    end
                    push!(candidates[root[ed]], cyc)
                    push!(used[nbr], tail)
                end
            end
        end
        push!(connectivity, Set(keys(pred)))
        # print(pred)
        setdiff!(nodeset, keys(pred))
    end
    cycles = []
    bicomps = []
    cyclemap = Dict(a.index => Set() for a in atomvector(mol))
    for bcmp in sort(collect(values(candidates)), by=length, rev=true)
        cycs = canonicalize_cycle.(minify_cycle(bcmp))
        cnt = length(cycles)
        for (i, c) in enumerate(cycs)
            for j in c
                push!(cyclemap[j], i + cnt)
            end
        end
        append!(cycles, cycs)
        push!(bicomps, collect(cnt + 1 : cnt + length(cycs)))
    end
    sort!(connectivity, by=length, rev=true)
    # TODO:cyclemap
    MolGraphTopology(cycles, bicomps, connectivity, cyclemap)
end

function molgraph_topology!(mol::MolecularGraph)
    mol.descriptor[:Topology] = molgraph_topology(mol)
    return
end


function canonicalize_cycle(nodes)
    """Re-order ring labels by index"""
    (fst, fstidx) = findmin(nodes)
    succidx = fstidx == lastindex(nodes) ? 1 : fstidx + 1
    succ = nodes[succidx]
    predidx = fstidx == 1 ? lastindex(nodes) : fstidx - 1
    pred = nodes[predidx]
    cp = succ < pred ? copy(nodes) : reverse(copy(nodes))
    while cp[1] != fst
        pushfirst!(cp, pop!(cp))
    end
    cp
end


function minify_cycle(bicomp; verbose=false)
    minified = []
    cnt = 0
    while !isempty(bicomp)
        cnt += 1
        if cnt > 100
            throw(OperationError("Cycle minimization failed"))
        end
        c = popfirst!(bicomp)
        init_c = c
        verbose ? print("$(length(c)) Cycle:$(c)\n") : nothing
        for m in minified
            verbose ? print("$(length(m)) Minified:$(m)\n") : nothing
            resolved = resolve_inclusion(c, m)
            if resolved != nothing
                if verbose
                    print(
                        "$(length(resolved[1])), $(length(resolved[2]))",
                        " Resolved:$(resolved[1]) $(resolved[2])\n"
                    )
                end
                c = resolved[1]
            end
        end
        verbose ? print("$(length(c)) New ring:$(c)\n") : nothing
        if length(c) == length(init_c) # no further minimization
            push!(minified, c)
        else
            push!(bicomp, c)
        end
    end
    minified
end


function resolve_inclusion(a, b)
    rev = length(a) > length(b)
    (lt, bg) = rev ? (copy(b), copy(a)) : (copy(a), copy(b))
    isec = intersect(lt, bg)
    isize = length(isec)
    if isize == 3 && length(lt) == 4
        # TODO: cubane special case
    elseif isize != length(lt) && isize <= length(lt) / 2 + 1
        return nothing # Already minimum
    end
    while !(bg[1] in isec) || bg[end] in isec
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
    new_bg = cat(bx, bg[1], reverse(lx), dims=1)
    rev ? (new_bg, lt) : (lt, new_bg)
end
