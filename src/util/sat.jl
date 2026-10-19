#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

# A small SAT solver (DPLL with unit propagation) for query containment.
# The formulas generated from query trees are small (typically less than 100 variables),
# so a simple pure Julia implementation is sufficient. This replaces PicoSAT.jl, whose
# binary (PicoSAT_jll) is not available on Windows.


"""
    is_satisfiable(clauses::Vector{Vector{Int}}) -> Bool

Return whether the CNF formula is satisfiable.

The formula is given in the DIMACS-like form: each clause is a vector of non-zero literals,
where a positive integer `i` is the variable `i` and a negative integer `-i` is its negation.
An empty clause is unsatisfiable and an empty formula is satisfiable.
"""
function is_satisfiable(clauses::Vector{Vector{Int}})
    nvars = 0
    for c in clauses, l in c
        nvars = max(nvars, abs(l))
    end
    assign = zeros(Int8, nvars)  # 1: true, -1: false, 0: unassigned
    return dpll!(clauses, assign)
end

# 1: true, -1: false, 0: unassigned
@inline function literal_value(assign::Vector{Int8}, l::Int)
    a = assign[abs(l)]
    return l > 0 ? a : -a
end

# Assign the literals of unit clauses until no unit clause remains.
# Return false if a conflict (a clause with all literals false) is found.
function unit_propagate!(clauses::Vector{Vector{Int}}, assign::Vector{Int8}, trail::Vector{Int})
    changed = true
    while changed
        changed = false
        for c in clauses
            unassigned = 0
            last = 0
            satisfied = false
            for l in c
                v = literal_value(assign, l)
                if v == 1
                    satisfied = true
                    break
                elseif v == 0
                    unassigned += 1
                    last = l
                end
            end
            satisfied && continue
            unassigned == 0 && return false  # conflict
            if unassigned == 1
                assign[abs(last)] = last > 0 ? Int8(1) : Int8(-1)
                push!(trail, abs(last))
                changed = true
            end
        end
    end
    return true
end

function dpll!(clauses::Vector{Vector{Int}}, assign::Vector{Int8})
    trail = Int[]
    if !unit_propagate!(clauses, assign, trail)
        for v in trail
            assign[v] = Int8(0)
        end
        return false
    end
    # branch on an unassigned variable of an unsatisfied clause
    var = 0
    for c in clauses
        any(l -> literal_value(assign, l) == 1, c) && continue
        for l in c
            if assign[abs(l)] == 0
                var = abs(l)
                break
            end
        end
        var == 0 || break
    end
    var == 0 && return true  # all clauses are satisfied
    for val in (Int8(1), Int8(-1))
        assign[var] = val
        dpll!(clauses, assign) && return true
    end
    assign[var] = Int8(0)
    for v in trail
        assign[v] = Int8(0)
    end
    return false
end
