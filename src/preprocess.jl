#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    kekulize!, kekulize,
    trivialhydrogens, allhydrogens,
    removehydrogens, addhydrogens,
    largestcomponentnodes, extractlargestcomponent,
    protonateacids!, protonateacids,
    deprotonateoniums!, deprotonateoniums,
    depolarize!, depolarize,
    polarize!, polarize,
    totriplebond!, totriplebond,
    toallenelike!, toallenelike


# TODO: large conjugated system
# TODO: salts and waters should detected by functional group analysis
# TODO: Phosphate, diphosphate, sulfate, nitrate, acetate,
# maleate, fumarate, succinate, citrate, tartrate, oxalate,
# mesylate, tosylate, besylate,
# benzoate, gluconate


"""
    kekulize!(mol::SMILES) -> Nothing

Kekulize the molecule that has SMILES aromatic bonds.

SMILES allows aromatic atoms in small letters - b, c, n, o, p, s, [as] and [se]. Once these are stored in `SmilesAtom.isaromatic` field, then `kekulize` will place double bonds to satisfy valences.

Kekulization is necessary for molecules parsed from SMILES. If not kekulized, some bond valence and implicit hydrogen properties would be wrong.
"""
function kekulize!(mol::SMILES)
    nodes = Set{Int}()
    for i in 1:nodecount(mol)
        nodeattr(mol, i).isaromatic === nothing && continue
        nodeattr(mol, i).isaromatic || continue
        if atomsymbol(mol)[i] === :C
            push!(nodes, i)
        elseif atomsymbol(mol)[i] in (:N, :P, :As)
            if (degree(mol, i) == 2 ||
                (degree(mol, i) == 3 && nodeattr(mol, i).charge == 1))
                push!(nodes, i)
            end
        end
    end
    subg = nodesubgraph(mol, nodes)
    coloring = twocoloring(subg)
    coloring === nothing && throw(
        ErrorException("Kekulization failed: Please check if your SMILES is valid (e.g. Pyrrole n should be [nH])"))
    a, b = coloring
    adjmap = Dict{Int,Set{Int}}(
        i => adjacencies(subg, i) for i in nodeset(subg))
    for (u, v) in maxcardmap(a, b, adjmap)
        e = findedgekey(mol, u, v)
        setedgeattr!(mol, e, setorder(edgeattr(mol, e), 2))
    end
    return mol
end


"""
    kekulize(mol::SMILES) -> SMILES

Return the kekulized molecule. See [`kekulize!`](@ref).
"""
function kekulize(mol::SMILES)
    newmol = graphmol(mol)
    kekulize!(newmol)
    return newmol
end



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

Return the molecule with hydrogen nodes removed.

If option `all` is set to true (default), all hydrogens will be removed, otherwise only trivial hydrogens will be removed (see [`trivialhydrogens`](@ref)).
"""
function removehydrogens(mol::GraphMol; all=true)
    hydrogens = all ? allhydrogens : trivialhydrogens
    ns = setdiff(nodeset(mol), hydrogens(mol))
    return graphmol(nodesubgraph(mol, ns))
end


"""
    addhydrogens(mol::GraphMol) -> GraphMol

Return the molecule with all hydrogen nodes explicitly attached.
"""
function addhydrogens(mol::GraphMol{A,B}) where {A<:Atom,B<:Bond}
    implicithconnected_ = implicithconnected(mol)
    newmol = graphmol(mol)
    for i in 1:nodecount(mol)
        for j in 1:implicithconnected_[i]
            addedge!(newmol, i, addnode!(newmol, A(:H)), B())
        end
    end
    return newmol
end



"""
    largestcomponentnodes(mol::GraphMol) -> Set{Int}

Return a set of nodes in the largest connected component.
"""
largestcomponentnodes(mol::GraphMol
    ) = sortstablemax(connectedcomponents(mol), by=length, init=Set{Int}())


"""
    extractlargestcomponent(mol::GraphMol) -> GraphMol

Return the largest connected component of the molecular graph.

This should be useful when you want to remove salt and water molecules from the molecular graph simply. On the other hand, this can remove important components from the mixture so carefully apply this preprocess method.
"""
extractlargestcomponent(mol::GraphMol
    ) = graphmol(nodesubgraph(mol, largestcomponentnodes(mol)))



"""
    protonateacids!(mol::GraphMol) -> GraphMol

Protonate oxo(thio) anion groups of the molecule.
"""
function protonateacids!(mol::GraphMol)
    atomsymbol_ = atomsymbol(mol)
    charge_ = charge(mol)
    connectivity_ = connectivity(mol)
    for o in findall(charge_ .== -1)
        atomsymbol_[o] in (:O, :S) || continue
        @assert connectivity_[o] == 1
        nbr = pop!(adjacencies(mol, o))
        charge_[nbr] == 1 && continue  # polarized double bond
        setnodeattr!(mol, o, setcharge(nodeattr(mol, o), 0))
    end
    clearcache!(mol)
end


"""
    protonateacids(mol::GraphMol) -> GraphMol

Return the molecule with its oxo(thio) anion groups protonated.

See [`protonateacids!`](@ref).
"""
function protonateacids(mol::GraphMol)
    newmol = graphmol(mol)
    protonateacids!(newmol)
    return newmol
end


"""
    deprotonateoniums!(mol::GraphMol) -> Nothing

Deprotonate onium groups of the molecule.
"""
function deprotonateoniums!(mol::GraphMol)
    hydrogenconnected_ = hydrogenconnected(mol)
    charge_ = charge(mol)
    for o in findall(charge_ .== 1)
        hydrogenconnected_[o] > 0 || continue
        setnodeattr!(mol, o, setcharge(nodeattr(mol, o), 0))
    end
    clearcache!(mol)
end


"""
    deprotonateoniums(mol::GraphMol) -> GraphMol

Return the molecule with its onium groups deprotonated.

See [`deprotonateoniums!`](@ref).
"""
function deprotonateoniums(mol::GraphMol)
    newmol = graphmol(mol)
    deprotonateoniums!(newmol)
    return newmol
end


"""
    depolarize!(mol::GraphMol; negative=:O, positive=[:C, :P]) -> Nothing

Depolarize dipole double bonds of the molecule.
"""
function depolarize!(mol::GraphMol; negative=:O, positive=[:C, :P])
    atomsymbol_ = atomsymbol(mol)
    charge_ = charge(mol)
    connectivity_ = connectivity(mol)
    isaromatic_ = isaromatic(mol)
    for o in findall(charge_ .== -1)
        atomsymbol_[o] === negative || continue
        @assert connectivity_[o] == 1
        (inc, adj) = iterate(neighbors(mol, o))[1]
        atomsymbol_[adj] in positive || continue
        charge_[adj] == 1 || continue
        isaromatic_[adj] && continue
        setnodeattr!(mol, o, setcharge(nodeattr(mol, o), 0))
        setnodeattr!(mol, adj, setcharge(nodeattr(mol, adj), 0))
        setedgeattr!(mol, inc, setorder(edgeattr(mol, inc), 2))
    end
    clearcache!(mol)
end


"""
    depolarize(mol::GraphMol) -> GraphMol

Return the molecule with its dipole double bonds depolarized.

See [`depolarize!`](@ref).
"""
function depolarize(mol::GraphMol; kwargs...)
    newmol = graphmol(mol)
    depolarize!(newmol; kwargs...)
    return newmol
end


"""
    polarize!(mol::GraphMol; negative=:O, positive=[:N, :S]) -> Nothing

Polarize dipole double bonds of the molecule.
"""
function polarize!(mol::GraphMol, negative=:O, positive=[:N, :S])
    atomsymbol_ = atomsymbol(mol)
    charge_ = charge(mol)
    bondorder_ = bondorder(mol)
    connectivity_ = connectivity(mol)
    for o in findall(connectivity_ .== 1)
        atomsymbol_[o] === negative || continue
        charge_[o] == 0 || continue
        (inc, adj) = iterate(neighbors(mol, o))[1]
        atomsymbol_[adj] in positive || continue
        charge_[adj] == 0 || continue
        bondorder_[inc] == 2 || continue
        setnodeattr!(mol, o, setcharge(nodeattr(mol, o), -1))
        setnodeattr!(mol, adj, setcharge(nodeattr(mol, adj), 1))
        setedgeattr!(mol, inc, setorder(edgeattr(mol, inc), 1))
    end
    clearcache!(mol)
end


"""
    polarize(mol::GraphMol) -> GraphMol

Return the molecule with its dipole double bonds polarized.

See [`polarize!`](@ref).
"""
function polarize(mol::GraphMol; kwargs...)
    newmol = graphmol(mol)
    polarize!(newmol; kwargs...)
    return newmol
end



function find13dipoles(mol::GraphMol)
    charge_ = charge(mol)
    pie_ = apparentvalence(mol) - nodedegree(mol)
    triads = Tuple{Int,Int,Int}[]
    for c in findall((pie_ .== 2) .* (charge_ .== 1))
        negs = Int[]
        notnegs = Int[]
        for (inc, adj) in neighbors(mol, c)
            push!(charge_[adj] == -1 ? negs : notnegs, adj)
        end
        (length(negs) == 1 && length(notnegs) == 1) || continue
        push!(triads, (negs[1], c, notnegs[1]))
    end
    return triads
end


"""
    totriplebond!(mol::GraphMol) -> Nothing

Standardize the molecule so that all 1,3-dipole groups are represented as triple bond and single bond (e.g. Diazo group C=[N+]=[N-] -> [C-][N+]#N).
"""
function totriplebond!(mol::GraphMol)
    for (first, center, third) in find13dipoles(mol)
        fe = findedgekey(mol, first, center)
        te = findedgekey(mol, center, third)
        edgeattr(mol, fe).order == 2 || continue
        setnodeattr!(mol, first, setcharge(nodeattr(mol, first), 0))
        setnodeattr!(mol, third, setcharge(nodeattr(mol, third), -1))
        setedgeattr!(mol, fe, setorder(edgeattr(mol, fe), 3))
        setedgeattr!(mol, te, setorder(edgeattr(mol, te), 1))
    end
    clearcache!(mol)
end


"""
    totriplebond(mol::GraphMol) -> GraphMol

Return the molecule standardized as the same way as [`totriplebond!`](@ref).
"""
function totriplebond(mol::GraphMol)
    newmol = graphmol(mol)
    totriplebond!(newmol)
    return newmol
end


"""
    toallenelike!(mol::GraphMol) -> Nothing

Standardize the molecule so that all 1,3-dipole groups are represented as allene-like structure (e.g. Diazo group [C-][N+]#N -> C=[N+]=[N-]).
"""
function toallenelike!(mol::GraphMol)
    for (first, center, third) in find13dipoles(mol)
        fe = findedgekey(mol, first, center)
        te = findedgekey(mol, center, third)
        edgeattr(mol, fe).order == 1 || continue
        setnodeattr!(mol, first, setcharge(nodeattr(mol, first), 0))
        setnodeattr!(mol, third, setcharge(nodeattr(mol, third), -1))
        setedgeattr!(mol, fe, setorder(edgeattr(mol, fe), 2))
        setedgeattr!(mol, te, setorder(edgeattr(mol, te), 2))
    end
    clearcache!(mol)
end


"""
    toallenelike(mol::GraphMol) -> GraphMol

Return the molecule standardized as the same way as [`toallenelike!`](@ref).
"""
function toallenelike(mol::GraphMol)
    newmol = graphmol(mol)
    toallenelike!(newmol)
    return newmol
end
