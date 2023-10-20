############################################################
# 
# Common utilities for rings of differential operators
# 
############################################################

coef_ring(D::AbstractDORing) = D |> unwrap |> base_ring

############################################################
# 
# Common utilities for differential operators
# 
############################################################

unwrap(wae::T) where T <: AbstractDiffOp = wae.elem

Base.parent(wae::T) where T <: AbstractDiffOp = wae.parent
# gens(wae::T) where T <: AbstractDiffOp = wae |> parent |> gens
# dgens(wae::T) where T <: AbstractDiffOp = wae |> parent |> dgens

function Base.show(io::IO, wae::T) where T <: AbstractDiffOp
    show(io, unwrap(wae))
end

Base.one(wae::Union{Type{T}, T}) where T <: AbstractDiffOp = T(parent(wae), one(unwrap(wae)))
Base.zero(wae::Union{Type{T}, T}) where T <: AbstractDiffOp = T(parent(wae), zero(unwrap(wae)))

function vars(wae::RatFuncElem)
    set = Set{typeof(wae)}()
    wae_nume = numerator(wae)
    wae_nume_vars = vars(wae_nume)
    set = union(set, wae_nume_vars)

    wae_deno = denominator(wae)
    wae_deno_vars = vars(wae_deno)
    set = union(set, wae_deno_vars)

    return set
end

"""
    vars(wae::T) where T <: AbstractDiffOp

# Examples 

```jldoctest
julia> D, (x, y), (dx, dy) = weyl_algebra(["x", "y"])
(2-d Weyl algebra in [x,y], PfaffianSystems.WAlgElem{AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.MPoly{Rational{BigInt}}}}[x, y], PfaffianSystems.WAlgElem{AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.MPoly{Rational{BigInt}}}}[dx, dy])
julia> vars(x+y)
2-element Vector{PfaffianSystems.WAlgElem{AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.MPoly{Rational{BigInt}}}}}:
 x
 y
```
"""
function vars(wae::T) where T<: AbstractDiffOp
    wae_coeffs = wae |> unwrap |> coefficients
    set = Set{typeof(wae)}()
    for c in wae_coeffs
        c_vars = vars(c)
        c_vars = unwrap(parent(wae)).(c_vars)
        c_vars = c_vars .|> (s->T(parent(wae), s))
        set = union(set, c_vars)

    end
    wae_vars = collect(set)
    
    re_wae_vars = Vector{typeof(wae)}()
    for i in gens(parent(wae))
        if i in wae_vars
            push!(re_wae_vars, i)
        end
    end

    return re_wae_vars
end


"""
    dvars(wae::T) where T <: AbstractDiffOp

# Examples

```jldoctest
julia> D, (x, y), (dx, dy) = weyl_algebra(["x", "y"])
(2-d Weyl algebra in [x,y], PfaffianSystems.WAlgElem{AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.MPoly{Rational{BigInt}}}}[x, y], PfaffianSystems.WAlgElem{AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.MPoly{Rational{BigInt}}}}[dx, dy])
julia> dvars(dx^2+y)
1-element Vector{PfaffianSystems.WAlgElem{AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.MPoly{Rational{BigInt}}}}}:
 dx
```
"""
function dvars(wae::T) where T <: AbstractDiffOp
    v = vars(unwrap(wae))
    wae_dvars = collect(Set(v .|> (s->T(parent(wae), s))))

    re_wae_dvars = Vector{typeof(wae)}()
    for i in dgens(parent(wae))
        if i in wae_dvars
            push!(re_wae_dvars, i)
        end
    end

    return re_wae_dvars 
end

isvar(dop::T) where T <: AbstractDiffOp = dop in gens(parent(dop))
isdvar(dop::T) where T <: AbstractDiffOp = dop in dgens(parent(dop))

"""

    evaluate(dop::T, vars::Vector{T}, vals::Vector{T}) where T <: AbstractDiffOp

Evaluate a differential operator with respect to a set of variables and their corresponding values.

# Arguments
- `dop::T`: A differential operator to evaluate.
- `vars::Vector{T}`: A vector of variables to evaluate the differential operator with respect to.
- `vals::Vector{T}`: A vector of values corresponding to the variables.

# Returns
- An instance of `T` representing the result of evaluating the differential operator.

The `vars` and `vals` vectors must have the same length, and each element in `vars` must be of the same type as its corresponding element in `vals`.

# Examples
```jldoctest
julia> D, (x, y), (dx, dy) = weyl_algebra(["x", "y"])
(2-d Weyl algebra in [x,y], PfaffianSystems.WAlgElem{AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.MPoly{Rational{BigInt}}}}[x, y], PfaffianSystems.WAlgElem{AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.MPoly{Rational{BigInt}}}}[dx, dy])
julia> evaluate(x*dx+y, [x], [x + y])
(x + y)*dx + y
julia> evaluate(x*dx + y, [dx], [dx+dy])
x*dx + x*dy + y
```
"""
function evaluate(dop::T, vrs::Vector{T}, vls::Vector{T}) where T <: AbstractDiffOp
    coefs = coefficients(dop) |> collect
    mons = monomials(dop) |> collect
    m_ring = parent(mons[1])

    v_pairs = [Vector{typeof(coefs[1])}(undef, 0), Vector{typeof(coefs[1])}(undef, 0)]
    dv_pairs = [Vector{typeof(mons[1])}(undef, 0), Vector{typeof(mons[1])}(undef, 0)]

    v2c(v) = collect(coefficients(unwrap(v)))[1] 
    v2m(v) = unwrap(v)


    for (var, val) in zip(vrs, vls)
        if isvar(var) &&  isempty(dvars(val))
            push!(v_pairs[1], v2c(var))
            push!(v_pairs[2], v2c(val))
        elseif isdvar(var) && isempty(vars(val)) 
            push!(dv_pairs[1], v2m(var))
            push!(dv_pairs[2], v2m(val))
        else
            if isvar(var)
                throw(DomainError("$var is variable while $val contains differential operator"))
            else
                throw(DomainError("$var is differential operator while $val contains variable"))
            end
        end
    end

    if length(v_pairs) > 0 
        coefs = [evaluate(c, v_pairs[1], v_pairs[2]) for c in coefs]
    end
    if length(dv_pairs) > 0 
        mons = [evaluate(m, dv_pairs[1], dv_pairs[2]) for m in mons]
    end

    return T(parent(dop), sum([m_ring(c)*m for (c, m) in zip(coefs, mons)]))
end

############################################################
# 
# Common arithmetic operations for differential operators
# 
############################################################

Base.:-(x::T) where T <: AbstractDiffOp = T(parent(x), -unwrap(x))

Base.:(==)(x::T, y::T) where T <: AbstractDiffOp  = unwrap(x) == unwrap(y)
Base.hash(x::T, h::UInt) where T <: AbstractDiffOp = hash(unwrap(x), h)

Base.:+(x::Union{Rational, Integer}, y::T) where T <: AbstractDiffOp  = T(parent(y), x + unwrap(y))
Base.:+(x::T, y::Union{Rational, Integer}) where T <: AbstractDiffOp = T(parent(x), unwrap(x) + y)
Base.:-(x::Union{Rational, Integer}, y::T) where T <: AbstractDiffOp = T(parent(y), x - unwrap(y))
Base.:-(x::T, y::Union{Rational, Integer}) where T <: AbstractDiffOp = T(parent(x), unwrap(x) - y)
Base.:*(x::Union{Rational, Integer}, y::T) where T <: AbstractDiffOp = T(parent(y), x * unwrap(y))
Base.:*(x::T, y::Union{Rational, Integer}) where T <: AbstractDiffOp = T(parent(x), unwrap(x) * y)

Base.:+(x::T, y::T) where T <: AbstractDiffOp = T(parent(x), unwrap(x) + unwrap(y))
Base.:-(x::T, y::T) where T <: AbstractDiffOp = T(parent(x), unwrap(x) - unwrap(y))

function Base.:*(l::T, r::T) where T <: AbstractDiffOp

    l == zero(parent(l)) && return zero(parent(l))
    r == zero(parent(r)) && return zero(parent(r))

    l_coeffs = l |> unwrap |> coefficients
    l_mons = l |> unwrap |> monomials
    r_coeffs = r |> unwrap |> coefficients
    r_mons = r |> unwrap |> monomials

    ret_dop = 0
    for (lc, lm) in zip(l_coeffs, l_mons), (rc, rm) in zip(r_coeffs, r_mons)
        ret_dop +=  lc * leibniz_rule(lm, rc) * rm
    end
    return T(parent(l), ret_dop)
end

function diff_op_pow(x::T, y::Integer) where T <: AbstractDiffOp
	y == 0 && return one(x)

    ret_dop = x
    for _ = 1:y-1
        ret_dop *= x
    end
    return ret_dop
end

Base.:literal_pow(::typeof(^), x::AbstractDiffOp,::Val{y}) where y = x^y


function leibniz_rule(l_mons::T, r_coeffs::U) where {U, T <: MPoly{<:U}}

    ret_dop = r_coeffs * l_mons
    n = nvars(parent(l_mons))

    for i=1:n, (m, c) in zip(monomials(ret_dop), coefficients(ret_dop))
        ret_dop += apply_diff(m,c,i)
    end

    return ret_dop
end

function apply_diff(l_mons::T ,r_coeffs::U, i::Integer) where {U, T <: MPoly{<:U}}

    ret_dop = 0
    k = 1
    coeff_var = gens(parent(r_coeffs))[i]
    mon_var = gens(parent(l_mons))[i]

    while true
        a = _nth_derivative(r_coeffs, coeff_var,k) * _nth_derivative(l_mons, mon_var,k) * parent(r_coeffs)(1//factorial(big(k)))
        a == 0 && break
        ret_dop += a
        k += 1
    end

    return ret_dop
end

function _nth_derivative(f::T, x::T ,n::Integer) where T

    n == 0 && return f 

    for _=1:n
        f = derivative(f,x)
    end

    return f

end

function derivative(f::RatFuncElem ,x::RatFuncElem)

    f_nume = numerator(f)
    f_deno = denominator(f)
    x_nume = numerator(x)

    nume = parent(x)((derivative(f_nume,x_nume)*f_deno - f_nume*derivative(f_deno,x_nume)))
    deno = parent(x)((f_deno^2))

    ret_dop = nume // deno
    return ret_dop
end

"""

    coefficients(f::T) where T <: AbstractDiffOp

Extract the coefficients of a differential operator.

Arguments:
- `f::T`: A differential operator to extract the coefficients from.

Returns:
- An array of instances of `T` representing the coefficients of the differential operator.
"""
coefficients(f::T) where T <: AbstractDiffOp = f |> unwrap |> coefficients |> collect

"""

    monomials(f::T) where T <: AbstractDiffOp

Extract the monomials of a differential operator.

Arguments:
- `f::T`: A differential operator to extract the monomials from.

Returns:
- An array of instances of `T` representing the monomials of the differential operator.
"""
function monomials(f::T) where T <: AbstractDiffOp

    f_mons = f |> unwrap |> monomials |> collect
    # f_mons = collect(f_mons .|> (s->T(parent(f), s)))

    return f_mons
end

"""

    exponent_vectors(f::T) where T <: AbstractDiffOp

Extract the exponent vectors of a differential operator.

Arguments:
- `f::T`: A differential operator to extract the exponent vectors from.

Returns:
- An array of arrays representing the exponent vectors of the differential operator.
"""
exponent_vectors(f::T) where T <: AbstractDiffOp = f |> unwrap |> exponent_vectors |> collect
        