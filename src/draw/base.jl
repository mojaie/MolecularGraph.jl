#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    Canvas,
    Color,
    DRAW_SETTING,
    draw2d_annot!,
    atomcolor,
    atomvisible,
    bondnotation,
    terminal_double_bond!,
    double_bond_along_ring!,
    boundary,
    chargesign,
    atommarkup,
    atomhtml


abstract type Canvas end


struct Color
    r::Int
    g::Int
    b::Int
end


const DRAW_SETTING = Dict(
    :display_terminal_carbon => false,
    :atomcolor => Dict(
        :H => Color(0, 0, 0),
        :B => Color(128, 0, 0),
        :C => Color(0, 0, 0),
        :N => Color(0, 0, 255),
        :O => Color(255, 0, 0),
        :F => Color(0, 255, 0),
        :Si => Color(128, 64, 192),
        :P => Color(192, 0, 192),
        :S => Color(192, 192, 0),
        :Cl => Color(64, 192, 64),
        :As => Color(128, 0, 128),
        :Se => Color(128, 128, 0),
        :Br => Color(0, 192, 0),
        :I => Color(0, 128, 0)
    ),
    :defaultatomcolor => Color(0, 192, 192)
)


function draw2d_annot!(mol::Molecule, setting)
    required_annotation(mol, :Topology)
    required_annotation(mol, :Valence)
    if mol.attribute[:sourcetype] == :sdfile
        mol.v[:Coords2D] = zeros(Float64, atomcount(mol), 2)
        for (i, a) in enumerate(mol.graph.nodes)
            mol.v[:Coords2D][i, :] = a.sdf_coords[1:2]
        end
    else
        mol.v[:Coords2D] = coords2d(mol)
    end
    mol.v[:AtomColor] = atomcolor(setting).(mol.v[:Symbol])
    mol.v[:AtomVisible] = atomvisible(
        setting[:display_terminal_carbon]).(mol.v[:Symbol], mol.v[:NumBonds])
    mol.v[:BondNotation] = bondnotation(mol)
    return
end

draw2d_annot!(mol::Molecule) = draw2d_annot!(mol, copy(DRAW_SETTING))


function atomcolor(setting)
    symbol -> get(setting[:atomcolor], symbol, setting[:defaultatomcolor])
end

function atomvisible(termc)
    (symbol, numb) -> numb == 0 || symbol != :C || (termc && numb == 1)
end


function bondnotation(mol::Molecule)
    if mol.attribute[:sourcetype] == :sdfile
        notation = [b.sdf_notation for b in mol.graph.edges]
    else
        notation = zeros(Int, bondcount(mol))
    end
    terminal_double_bond!(notation, mol.graph.adjacency,
                          mol.v[:NumBonds], mol.v[:Valence])
    double_bond_along_ring!(notation, mol.graph.adjacency,
                            mol.annotation[:Topology].cycles,
                            mol.v[:Coords2D], mol.v[:BondOrder])
    notation
end


function terminal_double_bond!(vec, adj, numb, valence)
    for nbrs in adj[(numb .== 1) .* (valence .== 2)]
        vec[collect(values(nbrs))[1]] = 2
    end
    return
end


function double_bond_along_ring!(vec, adj, rings, coords, order)
    for ring in sort(rings, by=length, rev=true)
        vtcs = [vec2d(coords[n, :]) for n in ring]
        cw = isclockwise(vtcs)
        if cw === nothing
            continue
        end
        directed = cw ? ring : reverse(ring)
        push!(directed, directed[1])
        for i in 1: length(directed) - 1
            u = directed[i]
            v = directed[i + 1]
            b = adj[u][v]
            if order[b] == 2 && u < v
                vec[b] = 1
            end
        end
    end
    return
end


function boundary(mol::Molecule, coords)
    (left, right) = extrema(coords[:, 1])
    (bottom, top) = extrema(coords[:, 2])
    width = right - left
    height = top - bottom
    dists = []
    # Size unit
    for bond in mol.graph.edges
        u = vec2d(coords[bond.u, :])
        v = vec2d(coords[bond.v, :])
        d = norm(v - u)
        if d > 0.0001  # Remove overlapped
            push!(dists, d)
        end
    end
    if isempty(dists)
        long = max(width, height)
        if long > 0.0001
            unit = long / sqrt(atomcount(mol))
        else
            unit = 1
        end
    else
        unit = median(dists) # Median bond length
    end

    (top, left, width, height, unit)
end


function chargesign(charge)
    if charge == 0
        return ""
    end
    sign = charge > 0 ? "+" : "–" # en dash, not hyphen-minus
    num = abs(charge)
    num > 1 ? string(num, sign) : sign
end


function atommarkup(symbol, charge, hcount, direction,
                    substart, subend, supstart, supend)
    if hcount == 1
        hs = "H"
    elseif hcount > 1
        hs = string("H", substart, hcount, subend)
    else
        hs = ""
    end
    chg = charge == 0 ? "" : string(supstart, chargesign(charge), supend)
    txt = direction == :left ? [chg, hs, symbol] : [symbol, hs, chg]
    string(txt...)
end


function atomhtml(symbol, charge, hcount, direction)
    atommarkup(symbol, charge, hcount, direction,
               "<sub>", "</sub>", "<sup>", "</sup>")
end
