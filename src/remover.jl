#
# This file is a part of graphmol.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    remove_H,
    remove_H!,
    removeall_H,
    removeall_H!,
    removewater,
    removewater!,
    removesalt,
    removesalt!


function remove_H!(mol::MutableMol)
    for (i, a) in mol.graph.nodes
        # TODO: check stereo (SMILES, SDFile)
        if (a.symbol == :H && a.charge == 0 && a.multiplicity == 1
                && a.mass === nothing)
            unlinkatom!(mol, i)
        end
    end
    return
end

function remove_H(mol; use_deepcopy=true)
    mol1 = use_deepcopy ? deepcopy(mol) : copy(mol)
    remove_H!(mol1)
    mol1
end


function removeall_H!(mol::MutableMol)
    for (i, a) in mol.graph.nodes
        if a.symbol == :H
            unlinkatom!(mol, i)
        end
    end
    return
end

function removeall_H(mol; use_deepcopy=true)
    mol1 = use_deepcopy ? deepcopy(mol) : copy(mol)
    removeall_H!(mol1)
    mol1
end


function removewater!(mol::MapMol)
    for (i, a) in mol.graph.nodes
        if a.symbol == :O && neighborcount(mol, i) == 0
            unlinkatom!(mol, i)
        end
    end
    return
end

function removewater(mol; use_deepcopy=true)
    mol1 = use_deepcopy ? deepcopy(mol) : copy(mol)
    removewater!(mol1)
    mol1
end


const SINGLE_ELEM_SALT = [:N, :Na, :Mg, :Al, :Cl, :K, :Ca, :Br, :I]


function removesalt!(mol::MapMol)
    for (i, a) in mol.graph.nodes
        if a.symbol in SINGLE_ELEM_SALT && neighborcount(mol, i) == 0
            unlinkatom!(mol, i)
        end
        # TODO: Phosphate, diphosphate, sulfate, nitrate, acetate,
        # maleate, fumarate, succinate, citrate, tartrate, oxalate,
        # mesylate, tosylate, besylate,
        # benzoate, gluconate
    end
    return
end

function removesalt(mol; use_deepcopy=true)
    mol1 = use_deepcopy ? deepcopy(mol) : copy(mol)
    removesalt!(mol1)
    mol1
end
