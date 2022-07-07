module PfaffianSystems

import Base: sort

# Write your package code here.
using Bijections
using Symbolics
using Symbolics: derivative, value
using SymbolicUtils: PolyForm
using DynamicPolynomials
using DynamicPolynomials: variables, exponents, coefficient, term
using DifferentialEquations: solve, ODEProblem

function Bijection{S, T}(dict::AbstractDict{S, T}) where S where T
	return Bijection(dict)
end

sort(v::Vector{Num}) = v[sortperm(string.(v))]
export sort

include("DiffOps.jl")
export genVars, genVars!, apply_do

using Symbolics: get_variables
using DataStructures: OrderedSet
export OrderedSet
include("DIdeals.jl")
export DIdeal, stdmon!, isZeroDimensional, makeTestVarsAndIdeal
export eliminationIdeal, intersectionIdeal, integrationIdeal

using Symbolics: scalarize
include("AsirWrapper.jl")
export isAsirAvailable, vec2str

include("PfaffSys.jl")
export PfaffianSystem, buildFuncA, integrate

end
