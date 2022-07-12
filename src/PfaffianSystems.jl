module PfaffianSystems

import Base: sort

# Write your package code here.
using Bijections
using Symbolics
# --- derivative of Symbolics does not keep coefficient type ---
# using Symbolics: derivative, value, get_variables, scalarize
# --------------------------------------------------------------
using Symbolics: value, get_variables, scalarize
using SymbolicUtils: PolyForm, BasicSymbolic
using DynamicPolynomials
using DynamicPolynomials: variables, exponents, coefficient, term
using DifferentialEquations: solve, ODEProblem
# using Symbolics: get_variables
using DataStructures: OrderedSet
# using Symbolics: scalarize
using Base: @invokelatest

function Bijection{S, T}(dict::AbstractDict{S, T}) where S where T
	return Bijection(dict)
end

sort(v::Vector{Num}) = v[sortperm(string.(v))]
export sort

include("AsirWrapper.jl")
export isAsirAvailable, vec2str, asir_derivative

include("DiffOps.jl")
export genVars, addVars, apply_do

export OrderedSet
include("DIdeals.jl")
export DIdeal, stdmon!, isZeroDimensional, makeTestVarsAndIdeal
export eliminationIdeal, intersectionIdeal, integrationIdeal, restrictionIdeal

include("PfaffSys.jl")
export PfaffianSystem, buildFuncA, integrate, applyStdMons

end
