module PfaffianSystems

import Base: sort

# Write your package code here.
using Bijections
using Symbolics
# --- derivative of Symbolics does not keep coefficient type ---
# using Symbolics: derivative, value, get_variables, scalarize
# --------------------------------------------------------------
using Symbolics: value, get_variables, scalarize, derivative
using SymbolicUtils: PolyForm, BasicSymbolic, isdiv, unpolyize
using DynamicPolynomials
using DynamicPolynomials: variables, exponents, coefficient, term
using DifferentialEquations: solve, ODEProblem
# using Symbolics: get_variables
using DataStructures: OrderedSet
# export OrderedSet
# using Symbolics: scalarize
using Base: @invokelatest

function Bijection{S, T}(dict::AbstractDict{S, T}) where S where T
	return Bijection(dict)
end

sort(v::Vector{Num}) = v[sortperm(string.(v))]
export sort

include("AsirWrapper.jl")
export isAsirAvailable, vec2str, asir_derivative, asir_reduce, asir_fctr

include("DiffOps.jl")
export genVars, addVars, apply_do

include("DIdeals.jl")
export DIdeal, stdmon!, isZeroDimensional, makeTestVarsAndIdeal, apply_ideal
export eliminationIdeal, intersectionIdeal, integrationIdeal, restrictionIdeal

include("PfaffSys.jl")
export PfaffianSystem, buildFuncA, integratePf, applyStdMons, denomLCM

end
