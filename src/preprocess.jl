#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    trivialhydrogens, allhydrogens,
    removehydrogens, addhydrogens!,
    largestcomponentnodes, extractlargestcomponent,
    neutralizeacids!, neutralizeoniums!, depolarize!, toallenelike!,
    kekulize!,
    preprocess!

# TODO: large conjugated system
# TODO: salts and waters should detected by functional group analysis
# TODO: Phosphate, diphosphate, sulfate, nitrate, acetate,
# maleate, fumarate, succinate, citrate, tartrate, oxalate,
# mesylate, tosylate, besylate,
# benzoate, gluconate


"""
    trivialhydrogens(mol::GraphMol) -> Set{Int}

Return a set of trivial hydrogen nodes.

"Trivial" means not charged, no unpaired electron, no specific mass, non-stereospecific and attached to organic heavy atoms. Note that stereochemstry of SDFile molecules will not be considered until `setstereocenter!` is applied.
"""
function trivialhydrogens(mol::GraphMol)
    hs = Set{Int}()
    organic_heavy = (
        :B, :C, :N, :O, :F, :Si, :P, :S, :Cl, :As, :Se, :Br, :I)
    for (i, atom) in enumerate(nodeattrs(mol))
        atom.symbol == :H || continue
        atom.charge == 0 || continue
        atom.multiplicity == 1 || continue
        atom.mass === nothing || continue
        @assert atom.stereo === :unspecified
        degree(mol, i) == 1 || continue
        (inc, adj) = iterate(neighbors(mol, i))[1]
        nodeattr(mol, adj).symbol in organic_heavy || continue
        edgeattr(mol, inc).order == 1 || continue
        nodeattr(mol, adj).stereo === :unspecified || continue
        if nodeattr(mol, adj) isa SmilesAtom
            @assert !nodeattr(mol, adj).isaromatic
        end
        push!(hs, i)
    end
    return hs
end


"""
    allhydrogens(mol::GraphMol) -> Set{Int}

Return a set of all hydrogen nodes.
"""
function allhydrogens(mol::GraphMol)
    hs = Set{Int}()
    for (i, atom) in enumerate(nodeattrs(mol))
        atom.symbol == :H || continue
        push!(hs, i)
    end
    return hs
end


"""
    removehydrogens(mol::GraphMol) -> GraphMol

A convenient method that returns the molecule with hydrogen nodes removed.

If option `all` is set to true (default), all hydrogens will be removed, otherwise only trivial hydrogens will be removed (see [`trivialhydrogens`](@ref)).

"""
function removehydrogens(mol::GraphMol; all=true)
    hydrogens = all ? allhydrogens : trivialhydrogens
    ns = setdiff(nodeset(mol), hydrogens(mol))
    return graphmol(nodesubgraph(mol, ns))
end


"""
    addhydrogens!(mol::GraphMol) -> GraphMol

Add hydrogen nodes to the molecular graph explicitly.

Note that this function edits `Atom` and `Bond` object fields directly (see [`Graph.clone`](@ref) and [`Graph.clearcache!`](@ref).
"""
function addhydrogens!(mol::GraphMol)
    implicithcount_ = implicithcount(mol)
    for i in 1:nodecount(mol)
        for j in 1:implicithcount_[i]
            n = addnode!(mol, nodeattrtype(mol)(:H))
            addedge!(mol, i, n, edgeattrtype(mol)())
        end
    end
end


"""
    largestcomponentnodes(mol::GraphMol) -> Set{Int}

Return a set of nodes in the largest connected component.
"""
largestcomponentnodes(mol::GraphMol
    ) = sortstablemax(connected_components(mol), by=length, init=Set{Int}())


"""
    extractlargestcomponent(mol::GraphMol) -> SubgraphView

Return largest connected component of the molecular graph.

This should be useful when you want to remove salt and water molecules from the molecular graph simply. On the other hand, this can remove important components from the mixture so carefully apply this preprocess method.
"""
extractlargestcomponent(mol::GraphMol
    ) = graphmol(nodesubgraph(mol, largestcomponentnodes(mol)))


"""
    neutralizeacids!(mol::GraphMol)

Neutralize oxo(thio) acids.

Note that this function edits `Atom` object fields directly (see [`Graph.clone`](@ref) and [`Graph.clearcache!`](@ref).
"""
function neutralizeacids!(mol::GraphMol)
    atomsymbol_ = atomsymbol(mol)
    charge_ = charge(mol)
    connectivity_ = connectivity(mol)
    pielectron_ = pielectron(mol)
    for o in findall(
            (atomsymbol_ .== :O) .* (charge_ .== -1) .* (connectivity_ .== 1))
        nbr = iterate(adjacencies(mol, o))[1]
        if pielectron_[nbr] == 1
            cnbrs = adjacencies(mol, nbr)
            pop!(cnbrs, o)
            for cn in cnbrs
                if (atomsymbol_[cn] in (:O, :S) && pielectron_[cn] == 1
                        && connectivity_[cn] == 1)
                    setnodeattr!(mol, o, setcharge(nodeattr(mol, o), 0))
                    break
                end
            end
        end
    end
end


"""
    neutralizeoniums!(mol::GraphMol)

Neutralize 1-3° oniums. Permanently charged quart-oniums will not be neutralized.

Note that this function edits `Atom` object fields directly (see [`Graph.clone`](@ref) and [`Graph.clearcache!`](@ref).
"""
function neutralizeoniums!(mol::GraphMol)
    for o in findall((charge(mol) .== 1) .* (hcount(mol) .> 0))
        setnodeattr!(mol, o, setcharge(nodeattr(mol, o), 0))
    end
end


"""
    depolarize!(mol::GraphMol)

Depolarize oxo groups except in the case that polarization is required for
aromaticity.

Note that this function edits `Atom` and `Bond` object fields directly (see [`Graph.clone`](@ref) and [`Graph.clearcache!`](@ref).
"""
function depolarize!(mol::GraphMol)
    charge_ = charge(mol)
    isaromatic_ = isaromatic(mol)
    for o in findall((atomsymbol(mol) .== :O) .* (charge_ .== -1))
        @assert degree(mol, o) == 1 "unexpected oxygen degree $(length(nbrs))"
        (inc, adj) = iterate(neighbors(mol, o))[1]
        if charge_[adj] == 1 && !isaromatic_[adj]
            setnodeattr!(mol, o, setcharge(nodeattr(mol, o), 0))
            setnodeattr!(mol, adj, setcharge(nodeattr(mol, adj), 0))
            setedgeattr!(mol, inc, setorder(edgeattr(mol, inc), 2))
        end
    end
end


"""
    toallenelike!(mol::GraphMol)

Convert triple bonds with charges into allene-like structure (ex. [C-][N+]#N -> C=[N+]=[N-]).

Note that this function edits `Atom` and `Bond` object fields directly (see [`Graph.clone`](@ref) and [`Graph.clearcache!`](@ref).
"""
function toallenelike!(mol::GraphMol)
    charge_ = charge(mol)
    for tb in findall(bondorder(mol) .== 3)
        tbond = edgeattr(mol, tb)
        (u, v) = getedge(mol, tb)
        for (f, s) in ((u, v), (v, u))
            nbrs = copy(neighbors(mol, f))
            pop!(nbrs, findedgekey(mol, f, s))
            length(nbrs) == 1 || continue
            (inc, adj) = iterate(nbrs)[1]
            if charge_[adj] == -1
                setnodeattr!(mol, adj, setcharge(nodeattr(mol, adj), 0))
                setnodeattr!(mol, s, setcharge(nodeattr(mol, s), -1))
                setedgeattr!(mol, inc, setorder(edgeattr(mol, inc), 2))
                setedgeattr!(mol, tb, setorder(edgeattr(mol, tb), 2))
            end
        end
    end
end


"""
    kekulize!(mol::SMILES)

Convert SMILES aromatic atoms (Organic atoms in lower case) into double bond representation.

Note that this function edits `Atom` and `Bond` object fields directly (see [`Graph.clone`](@ref) and [`Graph.clearcache!`](@ref).
"""
function kekulize!(mol::SMILES)
    nodes = Set{Int}()
    for i in 1:nodecount(mol)
        nodeattr(mol, i).isaromatic === true || continue
        if atomsymbol(mol)[i] === :C
            push!(nodes, i)
        elseif atomsymbol(mol)[i] in (:N, :P, :As) && hcount(mol)[i] == 0
            push!(nodes, i)
        end
    end
    subg = nodesubgraph(mol, nodes)
    a, b = twocoloring(subg)
    adjmap = Dict{Int,Set{Int}}(
        i => adjacencies(subg, i) for i in nodeset(subg))
    mapping = maxcardmap(a, b, adjmap)
    for (u, v) in mapping
        e = findedgekey(mol, u, v)
        setedgeattr!(mol, e, setorder(edgeattr(mol, e), 2))
    end
end


"""
    preprocess!(mol::GraphMol)

Default molecular preprocessing method.

- Ionized acids and oniums will be neutralized ([`neutralizeacids!`](@ref) and [`neutralizeoniums!`](@ref)).
- If there is several possible resonance structures, the one which have less total charge (if they are even, the one which have less maximum bond order) will be selected ([`depolarize!`](@ref) and [`toallenelike!`](@ref)).

For further preprocessing, following [`removehydrogens`](@ref), [`removestereohydrogens`](@ref) and [`largestcomponentgraph`](@ref) would work well.

Note that this function edits `Atom` and `Bond` object fields directly (see [`Graph.clone`](@ref) and [`Graph.clearcache!`](@ref)).
"""
function preprocess!(mol::GraphMol)
    neutralizeacids!(mol)
    neutralizeoniums!(mol)
    depolarize!(mol)
    toallenelike!(mol)
end
