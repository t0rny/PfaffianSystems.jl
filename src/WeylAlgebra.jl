########## References ##########
# GenericTypes.jl: https://github.com/Nemocas/AbstractAlgebra.jl/blob/597ce2fc0006a3ecdc24da2c63fc47ae0e98b7e6/src/generic/GenericTypes.jl
# 
################################

using AbstractAlgebra

const AA = AbstractAlgebra

############################################################
# 
#	WeylAlgebra / WAlgElem 
#
############################################################

# base_ring is the ring of coefficients
# s is the list of symbols
# ord is the ordering of the variables (not available yet, but will be)
# num_vars is the number of variables

# mutable struct WeylAlgebra{T <: RingElement} <: AA.NCRing
# 	base_ring::Ring
# 	s::Vector{Symbol}
# 	ord::Symbol
# 	num_vars::Int
# 	N::Int

# 	function WeylAlgebra{T}(R::Ring, s::Vector{Symbol}, ord::Symbol, cached::Bool = true) where T <: RingElement
# 		return AA.get_cached!(WAlgID, (R, s, ord), cached) do 
# 			new{T}(R, s, ord, length(s), 2*length(s))
# 		end::WeylAlgebra{T}
# 	end
# end

# function WeylAlgebra{T}(R::Ring, s::Vector{Symbol}, ord::Symbol, cached::Bool = true) where T <: RingElement
# 	@assert T == AA.elem_type(R)
# 	return WeylAlgebra{T}(R, s, ord, cached)
# end

# function WeylAlgebra{T}(R::Ring, s::Vector{Symbol}, cached::Bool = true) where T <: RingElement
# 	return WeylAlgebra{T}(R, s, :lex, cached)
# end

# const WAlgID = AA.CacheDictType{Tuple{Ring, Vector{Symbol}, Symbol}, Ring}()

# mutable struct WAlgElem{T <: RingElement} <: AA.NCRingElem
# 	coeffs::Vector{T}
# 	exps::Matrix{UInt}
# 	length::Int
# 	parent::WeylAlgebra{T}

# 	function WAlgElem{T}(D::WeylAlgebra) where T <: RingElement
# 		N = D.N
# 		return new{T}(Vector{T}(undef, 0), Matrix{UInt}(undef, N, 0), 0, D)
# 	end 

# 	WAlgElem{T}(D::WeylAlgebra, a::Vector{T}, b::Matrix{UInt}) where T <: RingElement = new{T}(a, b, length(a), D)

# 	function WAlgElem{T}(D::WeylAlgebra, a::T) where T <: RingElement
# 		N = D.N
# 		return iszero(a) ? new{T}(Vector{T}(undef, 0), Matrix{UInt}(undef, N, 0), 0, D) : new{T}([a], zeros(UInt, N, 1), 1, D)
# 	end
# end

WeylAlgebra{T} = MPolyRing{MPolyRingElem{T}} where T

############################################################
# 
# weyl_algebra constructor
# 
############################################################

"""
	WeylAlgebra
"""
# weyl_algebra(R::Ring, s::Union{Tuple{Vararg{T}}, AbstractVector{T}}; kw...)
function weyl_algebra(F::Field, s::Vector{Symbol}, ds::Vector{Symbol}; kw...)
	R, gens_R = polynomial_ring(F, s; kw...)
	D, gens_D = polynomial_ring(R, ds; kw...)
	(D, gens_R, gens_D)
end

function weyl_algebra(F::Field, s::AbstractVector{<:AbstractString}; kw...)
	weyl_algebra(F, Symbol.(s), Symbol.("d".*s); kw...)
end