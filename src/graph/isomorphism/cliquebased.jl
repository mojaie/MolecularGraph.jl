#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    nodemcsclique,
    edgemcsclique


function nodemcsclique(G, H;
        nodematcher=(g,h)->true, edgematcher=(g,h)->true, kwargs...)
    mcs = Dict{Int,Int}()
    for mapping in mcscliques(G, H, nodematcher, edgematcher; kwargs...)
        if length(mapping) > length(mcs)
            mcs = mapping
        end
    end
    return mcs
end


function edgemcsclique(G, H;
        nodematcher=(g,h)->true, edgematcher=(g,h)->true, kwargs...)
    lg = linegraph(G)
    lh = linegraph(H)
    nmatch = lgnodematcher(lg, lh, nodematcher, edgematcher)
    ematch = lgedgematcher(lg, lh, nodematcher)
    mcs = Dict{Int,Int}()
    for mapping in mcscliques(lg, lh, nmatch, ematch; kwargs...)
        if length(mapping) < length(mcs)
            continue  # cannot be a MCS
        end
        delta_y_correction!(mapping, G, H)
        if length(mapping) > length(mcs)
            mcs = mapping
        end
    end
    return mcs
end


function mcscliques(G, H, nodematcher, edgematcher; kwargs...)
    if nodecount(G) == 0 || nodecount(H) == 0
        return ()
    end
    if haskey(kwargs, :constraint)
        # TODO:
    else
        flt = modprodedgefilter(G, H, edgematcher)
    end
    prod = modularproduct(G, H, nodematcher, flt)
    clq = maximalcliques(prod; kwargs...)
    return map(clq) do nodes
        return Dict(nodeattr(prod, n).g => nodeattr(prod, n).h for n in nodes)
    end
end


function modprodedgefilter(G, H, edgematcher)
    return function (g1, g2, h1, h2)
        if hasedge(G, g1, g2) != hasedge(H, h1, h2)
            return false
        elseif !hasedge(G, g1, g2) && !hasedge(H, h1, h2)
            return true
        else
            return edgematcher(findedgekey(G, g1, g2), findedgekey(H, h1, h2))
        end
    end
end
