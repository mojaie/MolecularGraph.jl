#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    SMILESParser, SMARTSParser,
    smilestomol, smartstomol
    # resolvedefaultbond!


mutable struct SMILESParser{T,V,E}
    input::String
    pos::Int
    done::Bool
    node::Int # No. of current node
    branch::Int # No. of node at the current branch root
    root::Int # No. of node at the current tree root
    ringlabel::Dict{Int,Tuple{Int,Union{E,Symbol,Nothing}}} # TODO: strange union
    edges::Vector{Edge{T}}
    vprops::Vector{V}
    eprops::Vector{E}
end

SMILESParser{T}(smiles
    ) where T <: SimpleMolGraph = SMILESParser{eltype(T),vproptype(T),eproptype(T)}(
        smiles, 1, false, 0, 1, 1, Dict(), Edge{eltype(T)}[], vproptype(T)[], eproptype(T)[])



mutable struct SMARTSParser{T,V,E}
    input::String
    pos::Int
    done::Bool
    node::Int # No. of current node
    branch::Int # No. of node at the current branch root
    root::Int # No. of node at the current tree root
    ringlabel::Dict{Int,Tuple{Int,Union{E,Symbol,Nothing}}} # TODO: strange union
    edges::Vector{Edge{T}}
    vprops::Vector{V}
    eprops::Vector{E}
    connectivity::Vector{Vector{T}}
end

SMARTSParser{T}(smarts
    ) where T <: SimpleMolGraph = SMARTSParser{eltype(T),vproptype(T),eproptype(T)}(
        smarts, 1, false, 0, 1, 1, Dict(), Edge{eltype(T)}[], vproptype(T)[], eproptype(T)[], [])


"""
    smilestomol(smiles::AbstractString) -> GraphMol{SmilesAtom,SmilesBond}

Parse SMILES string into `GraphMol` object.

The syntax of SMILES in this library follows both Daylight SMILES and OpenSMILES.

# References

1. OpenSMILES Specification http://opensmiles.org/spec/open-smiles.html
1. Daylight Tutorials https://www.daylight.com/dayhtml_tutorials/index.html
"""
function smilestomol(::Type{T}, smiles::AbstractString) where T <: AbstractMolGraph
    state = SMILESParser{T}(smiles)
    fragment!(state)
    mol = T(state.edges, state.vprops, state.eprops)
    #kekulize && kekulize!(mol)
    return mol
end

smilestomol(smiles::AbstractString) = smilestomol(MolGraph{Int,SMILESAtom,SMILESBond}, smiles)


"""
    smartstomol(smarts::AbstractString) -> QueryMol{SmartsAtom,SmartsBond}

Parse SMARTS string into `QueryMol` object.
"""
function smartstomol(::Type{T}, smarts::AbstractString) where T <: AbstractMolGraph
    state = SMARTSParser{T}(smarts)
    occursin('.', smarts) ? componentquery!(state) : fragment!(state)
    mol = T(
        state.edges, state.vprops, state.eprops, Dict(:connectivity => state.connectivity))
    # setcache!(mol, :minimumcyclebasisnodes)
    # resolvedefaultbond!(mol)
    # mol = inferaromaticity(removehydrogens(mol))
    # kekulize && kekulize!(mol)
    # setdiastereo && setdiastereo!(mol)
    return mol
end

smartstomol(smarts::AbstractString
    ) = smartstomol(MolGraph{Int,QueryTree,QueryTree}, smarts)



# internal parser methods

function lookahead(state::Union{SMILESParser,SMARTSParser}, pos::Int)
    # Negative pos can be used
    newpos = state.pos + pos
    if newpos > length(state.input) || newpos < 1
        return Char(0)
    else
        c = state.input[newpos]
        if isascii(c)
            return state.input[newpos]
        else
            throw(ErrorException("invalid charactor $(c)"))
        end
    end
end

Base.read(state) = lookahead(state, 0)  # TODO: rename


function forward!(state::Union{SMILESParser,SMARTSParser}, num::Int)
    # Negative num can be used
    state.pos += num
    if state.pos > length(state.input)
        state.done = true
    elseif state.done
        throw(ErrorException("charactors in the buffer were used up"))
    elseif state.pos < 1
        throw(ErrorException("no more backtracking!"))
    end
end

forward!(state) = forward!(state, 1)
backtrack!(state, num) = forward!(state, -num)
backtrack!(state) = backtrack!(state, 1)



"""
    resolvedefaultbond(qmol::QueryMol) -> Nothing

Resolve default SMARTS bonds.
function resolvedefaultbond!(qmol::QueryMol)
    aromfml = QueryFormula(:isaromatic, true)
    singlefml = QueryFormula(:and, Set([
        QueryFormula(:order, 1),
        QueryFormula(:isaromatic, false)
    ]))
    notaromsymfml = QueryFormula(:and, Set([
        QueryFormula(:not, QueryFormula(:symbol, s)) for s in [:B, :C, :N, :O, :P, :S, :As, :Se]
    ]))
    # extract explicitly aromatic bonds
    arombonds = Set([])
    for ring in minimumcyclebasisnodes(qmol)
        isarom = true
        for n in ring
            nq = nodeattr(qmol, n).query
            if !issubset(nq, QueryFormula(:isaromatic, true), eval_recursive=false)
                isarom = false
                break
            end
        end
        if isarom
            union!(arombonds, edgeset(nodesubgraph(qmol, ring)))
        end
    end
    # CC, cC, Cc, [#6]C, C[#6], A*, *A -> -
    # cc, c[#6], [#6]c, [#6][#6] -> -,:
    # cc in the same aromatic ring -> :
    for e in 1:edgecount(qmol)
        q = edgeattr(qmol, e).query
        q == QueryFormula(:defaultbond, true) || continue
        if e in arombonds
            setedgeattr!(qmol, e, SmartsBond(aromfml))
            continue
        end
        (u, v) = getedge(qmol, e)
        uq = nodeattr(qmol, u).query
        unotarom = (issubset(uq,
            QueryFormula(:isaromatic, false), eval_recursive=false)
            || issubset(uq, notaromsymfml, eval_recursive=false))
        vq = nodeattr(qmol, v).query
        vnotarom = (issubset(vq,
            QueryFormula(:isaromatic, false), eval_recursive=false)
            || issubset(vq, notaromsymfml, eval_recursive=false))
        if unotarom || vnotarom
            setedgeattr!(qmol, e, SmartsBond(singlefml))
        else
            setedgeattr!(qmol, e, SmartsBond(QueryFormula(:or, Set([singlefml, aromfml]))))
        end
    end
end

"""