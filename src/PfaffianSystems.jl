module PfaffianSystems

# Write your package code here.
using Bijections
using Symbolics
using Symbolics: derivative, value
using SymbolicUtils: PolyForm
using DynamicPolynomials
using DynamicPolynomials: variables, exponents, coefficient, term


include("DiffOps.jl")
export genVars, genVars!, apply_do

using Symbolics: get_variables
using DataStructures: OrderedSet
export OrderedSet
include("DIdeals.jl")
export DIdeal, stdmon!, isZeroDimensional, makeTestVarsAndIdeal
export eliminationIdeal, idealIntersection

using Symbolics: scalarize
include("AsirWrapper.jl")
export isAsirAvailable, vec2str

include("PfaffSys.jl")
export PfaffianSystem, computePfaffSys

end
