#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    kekulize,
    trivialhydrogens, allhydrogens,
    removehydrogens, addhydrogens,
    largestcomponentnodes, extractlargestcomponent,
    protonateacids, deprotonateoniums,
    depolarize, polarize, totriplebond, toallenelike


# TODO: large conjugated system
# TODO: salts and waters should detected by functional group analysis
# TODO: Phosphate, diphosphate, sulfate, nitrate, acetate,
# maleate, fumarate, succinate, citrate, tartrate, oxalate,
# mesylate, tosylate, besylate,
# benzoate, gluconate


"""
    kekulize(mol::SMILES) -> SMILES

Return the molecule with SMILES aromatic bonds kekulized.

SMILES allows aromatic atoms in small letters - b, c, n, o, p, s, [as] and [se]. Once these are stored in `SmilesAtom.isaromatic` field, then `kekulize` will place double bonds to satisfy valences.

Kekulization should be applied for molecules parsed from SMILES. Otherwise, some bond valence and implicit hydrogen properties would be wrong.
"""
function kekulize(mol::SMILES)
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
    a, b = twocoloring(subg)
    adjmap = Dict{Int,Set{Int}}(
        i => adjacencies(subg, i) for i in nodeset(subg))
    mapping = maxcardmap(a, b, adjmap)
    newmol = graphmol(mol)
    for (u, v) in mapping
        e = findedgekey(mol, u, v)
        setedgeattr!(newmol, e, setorder(edgeattr(newmol, e), 2))
    end
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
    protonateacids(mol::GraphMol) -> GraphMol

Return the molecule with its oxo(thio) anion groups protonated.
"""
function protonateacids(mol::GraphMol)
    atomsymbol_ = atomsymbol(mol)
    charge_ = charge(mol)
    connectivity_ = connectivity(mol)
    newmol = graphmol(mol)
    for o in findall(charge_ .== -1)
        atomsymbol_[o] in (:O, :S) || continue
        @assert connectivity_[o] == 1
        nbr = pop!(adjacencies(mol, o))
        charge_[nbr] == 1 && continue  # polarized double bond
        setnodeattr!(newmol, o, setcharge(nodeattr(newmol, o), 0))
    end
    return newmol
end


"""
    deprotonateoniums(mol::GraphMol) -> GraphMol

Return the molecule with its onium groups are deprotonated.
"""
function deprotonateoniums(mol::GraphMol)
    hydrogenconnected_ = hydrogenconnected(mol)
    charge_ = charge(mol)
    newmol = graphmol(mol)
    for o in findall(charge_ .== 1)
        hydrogenconnected_[o] > 0 || continue
        setnodeattr!(newmol, o, setcharge(nodeattr(newmol, o), 0))
    end
    return newmol
end


"""
    depolarize(mol::GraphMol; negative=:O, positive=[:C, :P]) -> GraphMol

Return the molecule with its dipole double bonds are depolarized.
"""
function depolarize(mol::GraphMol; negative=:O, positive=[:C, :P])
    atomsymbol_ = atomsymbol(mol)
    charge_ = charge(mol)
    connectivity_ = connectivity(mol)
    isaromatic_ = isaromatic(mol)
    newmol = graphmol(mol)
    for o in findall(charge_ .== -1)
        atomsymbol_[o] === negative || continue
        @assert connectivity_[o] == 1
        (inc, adj) = iterate(neighbors(mol, o))[1]
        atomsymbol_[adj] in positive || continue
        charge_[adj] == 1 || continue
        isaromatic_[adj] && continue
        setnodeattr!(newmol, o, setcharge(nodeattr(newmol, o), 0))
        setnodeattr!(newmol, adj, setcharge(nodeattr(newmol, adj), 0))
        setedgeattr!(newmol, inc, setorder(edgeattr(newmol, inc), 2))
    end
    return newmol
end


"""
    polarize(mol::GraphMol; negative=:O, positive=[:N, :S]) -> GraphMol

Return the molecule with its dipole double bonds are polarized.
"""
function polarize(mol::GraphMol, negative=:O, positive=[:N, :S])
    atomsymbol_ = atomsymbol(mol)
    charge_ = charge(mol)
    bondorder_ = bondorder(mol)
    connectivity_ = connectivity(mol)
    newmol = graphmol(mol)
    for o in findall(connectivity_ .== 1)
        atomsymbol_[o] === negative || continue
        charge_[o] == 0 || continue
        (inc, adj) = iterate(neighbors(mol, o))[1]
        atomsymbol_[adj] in positive || continue
        charge_[adj] == 0 || continue
        bondorder_[inc] == 2 || continue
        setnodeattr!(newmol, o, setcharge(nodeattr(newmol, o), -1))
        setnodeattr!(newmol, adj, setcharge(nodeattr(newmol, adj), 1))
        setedgeattr!(newmol, inc, setorder(edgeattr(newmol, inc), 1))
    end
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
    totriplebond(mol::GraphMol) -> GraphMol

Return the molecule with its 1,3-dipole groups described in triple bond resonance structure (ex. Diazo group C=[N+]=[N-] -> [C-][N+]#N).
"""
function totriplebond(mol::GraphMol)
    newmol = graphmol(mol)
    for (first, center, third) in find13dipoles(mol)
        fe = findedgekey(mol, first, center)
        te = findedgekey(mol, center, third)
        edgeattr(mol, fe).order == 2 || continue
        setnodeattr!(newmol, first, setcharge(nodeattr(mol, first), 0))
        setnodeattr!(newmol, third, setcharge(nodeattr(mol, third), -1))
        setedgeattr!(newmol, fe, setorder(edgeattr(mol, fe), 3))
        setedgeattr!(newmol, te, setorder(edgeattr(mol, te), 1))
    end
    return newmol
end


"""
    toallenelike(mol::GraphMol) -> GraphMol

Return the molecule with its 1,3-dipole groups described in allene-like resonance structure (ex. Diazo group [C-][N+]#N -> C=[N+]=[N-]).
"""
function toallenelike(mol::GraphMol)
    newmol = graphmol(mol)
    for (first, center, third) in find13dipoles(mol)
        fe = findedgekey(mol, first, center)
        te = findedgekey(mol, center, third)
        edgeattr(mol, fe).order == 1 || continue
        setnodeattr!(newmol, first, setcharge(nodeattr(mol, first), 0))
        setnodeattr!(newmol, third, setcharge(nodeattr(mol, third), -1))
        setedgeattr!(newmol, fe, setorder(edgeattr(mol, fe), 2))
        setedgeattr!(newmol, te, setorder(edgeattr(mol, te), 2))
    end
    return newmol
end
