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

# WeylAlgebra{T} = MPolyRing{MPolyRingElem{T}} where T

struct WeylAlgebra{T <: MPolyRing{<:MPolyRingElem}}
	WAlg::T
end

function Base.show(io::IO, D::WeylAlgebra)
	print(io, "WAlg(", D.WAlg, ")")
end


struct WAlgElem{T <: MPolyRingElem}
	elem::T
end

function Base.show(io::IO, w::WAlgElem)
	print(io, w.elem)
end

Base.:+(x::WAlgElem, y::WAlgElem) = WAlgElem(x.elem + y.elem)


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
	(WeylAlgebra(D), WAlgElem.(D.(gens_R)), WAlgElem.(gens_D))
end

function weyl_algebra(F::Field, s::AbstractVector{<:AbstractString}; kw...)
	weyl_algebra(F, Symbol.(s), Symbol.("d".*s); kw...)
end