#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

export
    SmilesParser, SmartsParser,
    smilestomol, smartstomol,
    associate_operations, isequivalent, query_contains


mutable struct SmartsParserState{N<:AbstractNode,E<:UndirectedEdge}
    input::String
    allow_disconnected::Bool

    pos::Int
    done::Bool
    node::Int # No. of current node
    branch::Int # No. of node at the current branch root
    root::Int # No. of node at the current tree root
    ringlabel::Dict{Int,Tuple{Int,Union{E,Symbol,Nothing}}} # TODO: strange union

    edges::Vector{Tuple{Int,Int}}
    nodeattrs::Vector{N}
    edgeattrs::Vector{E}
    connectivity::Vector{Vector{Int}}

    function SmartsParserState{N,E}(str, disconn
            ) where {N<:AbstractNode,E<:UndirectedEdge}
        new(str, disconn, 1, false, 0, 1, 1, Dict(), [], [], [], [])
    end
end

SmilesParser = SmartsParserState{SmilesAtom,SmilesBond}
SmartsParser = SmartsParserState{SmartsAtom,SmartsBond}


function Base.parse(::Type{SMILES}, str::AbstractString)
    state = SmilesParser(str, true)
    fragment!(state)
    return graphmol(state.edges, state.nodeattrs, state.edgeattrs)
end

function Base.parse(::Type{SMARTS}, str::AbstractString)
    if occursin('.', str)
        state = SmartsParser(str, true)
        componentquery!(state)
    else
        state = SmartsParser(str, false)
        fragment!(state)
    end
    qmol = querymol(state.edges, state.nodeattrs, state.edgeattrs, state.connectivity)
    resolvedefaultbond!(qmol)
    return qmol
end


"""
    smilestomol(smiles::AbstractString) -> GraphMol{SmilesAtom,SmilesBond}

Parse SMILES string into `GraphMol` object.
"""
function smilestomol(smiles::AbstractString)
    mol = parse(SMILES, smiles)
    kekulize!(mol)
    setdiastereo!(mol)
    return mol
end


"""
    smartstomol(smarts::AbstractString) -> QueryMol{SmartsAtom,SmartsBond}

Parse SMARTS string into `QueryMol` object.
"""
smartstomol(smarts::AbstractString) = parse(SMARTS, smarts)



function lookahead(state::SmartsParserState, pos::Int)
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


function forward!(state::SmartsParserState, num::Int)
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
    associate_operations(fml::Pair) -> Pair

Return juxtaposed formulae (ex. :and => (A, :and => (B, C)) -> :and => (A, B, C)).
"""
function associate_operations(fml::Pair)
    fml.first in (:and, :or) || return fml
    associated = Set()
    for elem in fml.second
        if elem.first === fml.first
            union!(associated, associate_operations(elem).second)
        else
            push!(associated, elem)
        end
    end
    return fml.first => Tuple(associated)
end


"""
    isequivalent(fml1::Pair, fml2::Pair) -> Pair

Check if the two formulae are equivalent.
"""
function isequivalent(fml1::Pair, fml2::Pair)
    fml1.first === fml2.first || return false
    fml1.first in (:and, :or, :not) || return fml1.second == fml2.second
    fml1.first === :not && return isequivalent(fml1.second, fml2.second)
    length(fml1.second) == length(fml2.second) || return false
    f1map = Dict(i => v for (i, v) in enumerate(fml1.second))
    f2map = Dict(i => v for (i, v) in enumerate(fml2.second))
    iseq = (x, y) -> isequivalent(f1map[x], f2map[y])
    return maxcard(keys(f1map), keys(f2map), iseq) == length(f1map)
end


"""
    query_contains(fml1::Pair, fml2::Pair) -> Pair

Check if the two formulae are equivalent.
"""
function query_contains(fml1::Pair, fml2::Pair)
    fml1 == (:any => true) && return true
    if fml2.first === :and
        f2map = Dict(i => v for (i, v) in enumerate(fml2.second))
        if fml1.first === :and
            f1map = Dict(i => v for (i, v) in enumerate(fml1.second))
        else
            f1map = Dict(1 => fml1)
        end
        func = (x, y) -> query_contains(f1map[x], f2map[y])
        return maxcard(keys(f1map), keys(f2map), func) == length(f1map)
    elseif fml1.first === :or
        f1map = Dict(i => v for (i, v) in enumerate(fml1.second))
        if fml2.first === :or
            f2map = Dict(i => v for (i, v) in enumerate(fml2.second))
        else
            f2map = Dict(1 => fml2)
        end
        func = (x, y) -> query_contains(f1map[x], f2map[y])
        return maxcard(keys(f1map), keys(f2map), func) == length(f2map)
    else
        return fml1 == fml2
    end
end


"""
    resolvedefaultbond(qmol::QueryMol) -> Nothing

Resolve default SMARTS bonds.
"""
function resolvedefaultbond!(qmol::QueryMol)
    newedges = QueryBond[]
    for e in 1:edgecount(qmol)
        (u, v) = getedge(qmol, e)
        q = edgeattr(qmol, e).query
        if q == (:defaultbond => true)
            uq = nodeattr(qmol, u).query
            isuqalip = uq == (:isaromatic => false) || (uq.first === :and && (:isaromatic => false) in uq.second)
            isuqarom = uq == (:isaromatic => true) || (uq.first === :and && (:isaromatic => true) in uq.second)
            vq = nodeattr(qmol, v).query
            isvqalip = vq == (:isaromatic => false) || (uq.first === :and && (:isaromatic => false) in uq.second)
            isvqarom = vq == (:isaromatic => true) || (vq.first === :and && (:isaromatic => true) in vq.second)
            if isuqarom && isvqarom
                push!(newedges, SmartsBond(:isaromaticbond => true))
            elseif (isuqalip || isuqarom) && (isvqalip || isvqarom)
                push!(newedges, SmartsBond(:bondorder => 1))
            else
                push!(newedges, SmartsBond(:or => (:bondorder => 1, :isaromaticbond => true)))
            end
        else
            push!(newedges, edgeattr(qmol, e))
        end
    end
    empty!(qmol.edgeattrs)
    append!(qmol.edgeattrs, newedges)
end