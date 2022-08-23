module PfaffianSystems

import Base: sort

# Write your package code here.
using Base: @invokelatest
using Bijections
using Symbolics
# --- derivative of Symbolics does not keep coefficient type ---
# using Symbolics: derivative, value, get_variables, scalarize
# --------------------------------------------------------------
using Symbolics: value, get_variables, scalarize, derivative, wrap, unwrap
using SymbolicUtils: PolyForm, BasicSymbolic, isdiv, unpolyize
using DynamicPolynomials
using DynamicPolynomials: variables, exponents, coefficient, term
using DifferentialEquations: solve, ODEProblem
# using Symbolics: get_variables
using DataStructures: OrderedSet, OrderedDict
# export OrderedSet
# using Symbolics: scalarize
using MultivariatePolynomials

const DP = DynamicPolynomials
const MP = MultivariatePolynomials
const SU = SymbolicUtils


function Bijection{S, T}(dict::AbstractDict{S, T}) where {S, T}
	return Bijection(dict)
end
function empty!(b::Bijection{S, T}) where {S, T}
	dom = domain(b)
	for x in dom
		delete!(b, x)
	end
	return b
end

sort(v::Vector{Num}) = v[sortperm(string.(v))]
export sort

include("AsirWrapper.jl")
export isAsirAvailable, vec2str, asir_derivative, asir_reduce, asir_fctr

include("DiffOps.jl")
export PolyDiffOp
export genVars, addVars, apply_do, dmul, genVars2
# export makeTestEnvs, canonicalize, RatDiffOp, DGen

include("DIdeals.jl")
export DIdeal, stdmon!, isZeroDimensional, makeTestVarsAndIdeal, apply_ideal
export eliminationIdeal, intersectionIdeal, integrationIdeal, restrictionIdeal
export DIdeal2, makeTestVarsAndIdeal2

include("PfaffSys.jl")
export PfaffianSystem, get_vars, get_dvars, buildFuncA, integratePf, applyStdMons, denomLCM

end
