
############################################################
# 
#	WeylAlgebra / WAlgElem 
#
############################################################

struct WeylAlgebra{T <: MPolyRing{<:MPoly}} <: AbstractDORing
	WAlg::T

    function WeylAlgebra(F::Field, s::Vector{Symbol}; kw...) 
        length(s) != length(unique(s)) && throw(ArgumentError("variables must be unique"))
        ds = Symbol.("d" .* string.(s))
        # polyR = AbstractAlgebra.PolynomialRing_only(F, s; kw...)
        # raw_D = AbstractAlgebra.PolynomialRing_only(polyR, ds; kw...)
        polyR, _ = PolynomialRing(F, s; kw...)
        raw_D, _ = PolynomialRing(polyR, ds; kw...)
        new{typeof(raw_D)}(raw_D)
    end
end
unwrap(D::WeylAlgebra) = D.WAlg
(D::WeylAlgebra)(num::Union{Rational, Integer}) = WAlgElem(D, unwrap(D)(num))
function (D::WeylAlgebra)(poly::MPoly)
    !(parent(poly) == unwrap(D)) && throw(DomainError("$poly is not an element of $(unwrap(D))"))
    return WAlgElem(D, poly)
end

elem_type(D::Union{Type{WeylAlgebra{T}}, WeylAlgebra{T}}) where {S <: MPoly, T <: MPolyRing{S}} = WAlgElem{Generic.MPoly{S}}

function Base.show(io::IO, D::WeylAlgebra)
    print(io, nvars(D),"-d Weyl algebra in [$(join(string.(gens(D)), ","))]")
end

# TODO: wrapper for constant

struct WAlgElem{T <: MPoly{<:MPoly}} <: AbstractDiffOp
    parent::WeylAlgebra
	elem::T
end

function Base.:^(x::WAlgElem, y::Integer)
    y < 0 && return throw(DomainError("exponent must be non-negative"))
    return diff_op_pow(x,y)
end

############################################################
# 
# coersions
# 
############################################################

function genNo2str_dict(D::AbstractDORing)
    bj = Bijection{Integer, String}()
    for (i, g) in enumerate(gens(D))
        bj[i] = string(g)
    end
    return bj
end

function coercion_homomorphism(D1::AbstractDORing, D2::AbstractDORing)
    bj1 = genNo2str_dict(D1) 
    bj2 = genNo2str_dict(D2)

    common_vars_str = intersect(bj1.range, bj2.range)
    common_vars_str |> isempty && throw(DomainError("Cannot determine homomorphism from " * string(D1) * " to " * string(D2)))

    hom = Bijection{Integer, Integer}()

    for s in common_vars_str
        hom[bj1(s)] = bj2(s)
    end

    return hom
end



function _coerce_unsafe(x::MPoly, M::MPolyRing, index_map::Dict{Integer, <:Integer})
    n = length(index_map)
    cezip = zip(coefficients(x), exponent_vectors(x))
    C = Generic.MPolyBuildCtx(M)

    for (c,e) in cezip
        push_term!(C, c,[get(e, index_map[i], 0) for i in 1:n] )
    end

    return finish(C)
end

function coerce(x::WAlgElem, D::WeylAlgebra)
    hom = coercion_homomorphism(parent(x), D) 
    n = nvars(D)

    # define mapping from index of D to index of parent(x)
    index_map = Dict{Integer, Integer}()
    invhom = inv(hom)

    for i = 1:n
        index_map[i] = get(invhom, i, 0)
    end

    cezip = zip(coefficients(unwrap(x)), exponent_vectors(unwrap(x)))
    C = Generic.MPolyBuildCtx(unwrap(D))
    M = base_ring(D)

    for (c, e) in cezip
        coerced_c = _coerce_unsafe(c, M, index_map)
        push_term!(C, coerced_c, [get(e, index_map[i], 0) for i in 1:n])
    end

    WAlgElem(D, finish(C))
end

(D::WeylAlgebra)(x::WAlgElem) = coerce(x, D)

############################################################
# 
# Weyl algebra constructors
# 
############################################################

"""
	WeylAlgebra
"""
function weyl_algebra(F::Field, s::AbstractVector{<:AbstractString}; kw...)
    D = WeylAlgebra(F, Symbol.(s))
    return D, gens(D), dgens(D)
end
weyl_algebra(s::AbstractVector{<:AbstractString}; kw...) = weyl_algebra(QQ, s; kw...)

function weyl_algebra(F::Field, s::AbstractString; kw...)
    D = WeylAlgebra(F, [Symbol(s)])
    return D, gens(D)[1], dgens(D)[1]
end
weyl_algebra(s::AbstractString; kw...) = weyl_algebra(QQ, s; kw...)

function weyl_algebra(F::Field, D::WeylAlgebra, new_vars::AbstractVector{<:AbstractString}; kw...)
    old_vars = gens(D) .|> string
    if !(:new_vars_pos in keys(kw)) || kw[:new_vars_pos] == :append
        vars = [old_vars; new_vars]
    elseif kw[:new_vars_pos] == :prepend
        vars = [new_vars; old_vars]
    else
        throw(ArgumentError("new_vars_pos must be :prepend or :append"))
    end
    De, g, dg = weyl_algebra(F, vars; kw...)
    return De, g, dg
end
weyl_algebra(D::WeylAlgebra, s::AbstractVector{<:AbstractString}; kw...) = weyl_algebra(QQ, D, s; kw...)

function weyl_algebra(F::Field, s::AbstractString, n::Integer)
    D = WeylAlgebra(F, [Symbol(s,i) for i = 1:n])
    return D, gens(D), dgens(D)
end
weyl_algebra(s::AbstractString,n::Integer) = weyl_algebra(QQ, s, n)

# TODO: make new Weyl algebra with a part of variables


############################################################
# 
# Common to WeylAlgebra and DiffOpRing
# 
############################################################

Base.one(D::T) where T <: AbstractDORing = D(1)
Base.zero(D::T) where T <: AbstractDORing = D(0)

base_ring(D::AbstractDORing) = D |> unwrap |> base_ring

function gens(D::AbstractDORing)
    g = D |> base_ring |> gens
    g = unwrap(D).(g)
    # g .|> WAlgElem
    return g .|> (s->elem_type(D)(D, s))
end

function dgens(D::AbstractDORing)
    dg = D |> unwrap |> gens
    return dg .|> (s->elem_type(D)(D, s))    
end

nvars(D::AbstractDORing) = D |> unwrap |> nvars

Base.:(==)(x::AbstractDORing, y::AbstractDORing) = unwrap(x) == unwrap(y)