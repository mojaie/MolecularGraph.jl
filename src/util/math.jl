#
# This file is a part of MolecularGraph.jl
# Licensed under the MIT License http://opensource.org/licenses/MIT
#

"""
    logfactorial(n::Int) -> Float64

Return the approximate value of ln ``n!`` by Ramanujan method.
"""
function logfactorial(n::Int)
    n >= 0 || throw(DomainError(n, "n should be 0 or more"))
    n == 0 && return 0.0
    return n * log(n) - n + 1/6 * log(8*n^3 + 4*n^2 + n + 1/30) + 1/2 * log(pi)
end
