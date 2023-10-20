
############################################################
# 
# DiffOpRing / DiffOpRingElem
#	Ring of differential operators over rational functions
#
############################################################

struct DiffOpRing{T <: MPolyRing{<:RatFuncElem}} <: AbstractDORing
	DOR::T

    function DiffOpRing(F::Field, s::Vector{Symbol}; kw...)
        length(s) != length(unique(s)) && throw(ArgumentError("variables must be unique"))
        ds = Symbol.("d" .* string.(s))
        R, _ = RationalFunctionField(QQ, string.(s))
        raw_D = AbstractAlgebra.polynomial_ring_only(R, ds; kw...)
        new{typeof(raw_D)}(raw_D)
    end

end
unwrap(R::DiffOpRing) = R.DOR
(R::DiffOpRing)(num::Union{Rational, Integer}) = DORElem(R, unwrap(R)(num))
function (R::DiffOpRing)(poly::MPoly)
    !(parent(poly) == unwrap(R)) && throw(DomainError("$poly is not an element of $(unwrap(R))"))
    return DORElem(R, poly)
end


elem_type(D::Union{Type{DiffOpRing{T}}, DiffOpRing{T}}) where {S <: RatFuncElem, T <: MPolyRing{S}} = DORElem{Generic.MPoly{S}}

function Base.show(io::IO, R::DiffOpRing)
	print(io, nvars(R), "-dimensional ring of differential opeartors in [$(join(string.(gens(R)), ","))]")
end

struct DORElem{T <: MPoly{<:RatFuncElem}} <: AbstractDiffOp
    parent::DiffOpRing
	elem::T
end

function Base.:^(x::DORElem, y::Integer) 
    if y < 0
        x = 1 // x
        y = - y
    end
    return diff_op_pow(x,y)
end

Base.://(x::DORElem,y::Union{Rational, Integer}) = x//(parent(x)(y))
Base.://(x::Union{Rational, Integer},y::DORElem) = (parent(y)(x))//y

function Base.://(x::DORElem,y::DORElem)
    x_coeff = x |> unwrap |> coefficients                    
    x_mon = x |> unwrap |> monomials  
    y_coeff = y |> unwrap |> coefficients    
    y_mon = y |> unwrap |> monomials

    ret_dop = 0

    isempty(dvars(y)) != true && return throw(DomainError("division by differential operator", string(y)))

    for (xc, xm) in zip(x_coeff, x_mon), (yc, ym) in zip(y_coeff, y_mon)
        ret_dop += (xc // yc) * xm
    end

    return DORElem(parent(x), ret_dop)

end

############################################################
# 
# coersions
# 
############################################################

function _coerce_unsafe(x::RatFuncElem, R::Generic.RationalFunctionField, index_map::Dict{<:Integer, <:Integer})
    n = length(index_map)
    MPoly = R |> zero |> numerator |> parent
    #coercion of numerator
    x_nume = numerator(x)
    coerced_nume = _coerce_unsafe(x_nume, MPoly, index_map)
    
    #coercion of denominator
    x_deno = denominator(x)
    coerced_deno = _coerce_unsafe(x_deno, MPoly, index_map)
    
    # return coerced numerator divided ny coerced denominator

    return R(coerced_nume) // R(coerced_deno)
end

function _coerce_unsafe(x::MPoly, R::Generic.RationalFunctionField, index_map::Dict{<:Integer, <:Integer})
    n = length(index_map)
    MPoly = R |> zero |> numerator |> parent
    cezip = zip(coefficients(x), exponent_vectors(x))
    C = Generic.MPolyBuildCtx(MPoly)
    for (c,e) in cezip
        push_term!(C, c,[get(e, index_map[i], 0) for i in 1:n] )
    end
    return R(finish(C))
end


function coerce(x::DORElem, D::DiffOpRing)
    hom = coercion_homomorphism(parent(x), D)
    n = nvars(D)

    index_map = Dict{Integer, Integer}()
    invhom = inv(hom)
    for i = 1:n
        index_map[i] = get(invhom, i, 0)
    end

    cezip = zip(coefficients(unwrap(x)), exponent_vectors(unwrap(x)))
    M = Generic.MPolyBuildCtx(unwrap(D))

    R = base_ring(D)
    for (c,e) in cezip
        coerced_c = _coerce_unsafe(c, R ,index_map)
        push_term!(M, coerced_c, [get(e, index_map[i], 0) for i in 1:n])
    end
    DORElem(D, finish(M)) 
end
(R::DiffOpRing)(x::DORElem) = coerce(x, R)



function coerce(x::WAlgElem, D::DiffOpRing)
    hom = coercion_homomorphism(parent(x), D) 
    n = nvars(D)

    # define mapping from index of D to index of parent(x)
    index_map = Dict{Integer, Integer}()
    invhom = inv(hom)

    for i = 1:n
        index_map[i] = get(invhom, i, 0)
    end

    cezip = zip(coefficients(unwrap(x)), exponent_vectors(unwrap(x)))
    M = Generic.MPolyBuildCtx(unwrap(D))
    R = base_ring(D)

    for (c, e) in cezip
        coerced_c = _coerce_unsafe(c, R, index_map)
        push_term!(M, coerced_c, [get(e, index_map[i], 0) for i in 1:n])
    end

    DORElem(D, finish(M))
end
(R::DiffOpRing)(x::WAlgElem) = coerce(x, R)


############################################################
# 
# DiffOpRing constructor
# 
############################################################

"""
	Ring of differential operators over rational functions
"""
function diff_op_ring(F::Field, s::AbstractVector{<:AbstractString}; kw...)
	D = DiffOpRing(F, Symbol.(s))
    return D, gens(D), dgens(D)
end
diff_op_ring(s::AbstractVector{<:AbstractString}; kw...) = diff_op_ring(QQ, s; kw...)

function diff_op_ring(F::Field, s::AbstractString; kw...)
    D = DiffOpRing(F, [Symbol(s)])
    return D, gens(D)[1], dgens(D)[1]
end
diff_op_ring(s::AbstractString; kw...) = diff_op_ring(QQ, s; kw...)

function diff_op_ring(F::Field, s::AbstractString, n::Integer)
    D = DiffOpRing(F, [Symbol(s,i) for i = 1:n])
    return D, gens(D), dgens(D)
end
diff_op_ring(s::AbstractString,n::Integer) = diff_op_ring(QQ, s, n)


############################################################
# 
# Ideal of ring of differential operators
# 
############################################################

function grevlex_tie_break(x::Vector{<:Int}, y::Vector{<:Int})
    last_nz_idx = findlast(x-y .!= 0)
    (x-y)[last_nz_idx] < 0 && return true
    return false
end

"""
    leading_term(f::DORElem, order::Symbol=:lex)
Return the leading term of `f` with respect to `order`. 
Only `order=:lex` and `order=:grevlex` are supported now.
Default value of `order` is `:grevlex`.
"""
function leading_term(f::DORElem; order::Symbol=:grevlex)
    f == zero(parent(f)) && return zero(parent(f))
    f_coes = coefficients(f)
    f_mons = monomials(f)
    if order == :lex
        # f_mon = f_mons[1] |> unwrap
        return DORElem(parent(f), f_coes[1] * f_mons[1])
    elseif order == :revlex
    elseif order == :grlex
    elseif order == :grevlex
        exps = exponent_vectors(f)
        exp_sum = sum.(exps)
        exp_max = findall(x->x==maximum(exp_sum), exp_sum)
        sort!(exp_max, rev=true, lt = (x,y) -> grevlex_tie_break(exps[x], exps[y]))
        return DORElem(parent(f), f_coes[exp_max[1]] * f_mons[exp_max[1]])
        # for i in 1:size(dgens(parent(f)))[1]
        #     # for j in exp_max
        #     #     # exponent_vectors(f)[j][]


        #     # end
        # end
    end
end

"""
    normalform(f::T, G::Vector{T}) where T<: DORElem
Compute normal form of `f` with respect to `G` and return the remainder `r` and quotients `q`.
"""
function normalform(f::T, G::Vector{T}; order::Symbol=:grevlex) where T<: DORElem
    r = zero(parent(f))
    q = [zero(parent(f)) for _ in G]

    while f != zero(parent(f))
        r_1, q_1 = wnormalform(f, G, order=order)
        r = r + leading_term(r_1, order=order)

        for i in 1:size(q)[1]
            q[i] = q[i] + q_1[i]
        end

        f = r_1 - leading_term(r_1, order=order)
    end
    return r, q
end

"""
  wnormalform(f::T, G::Vector{T}) where T<: DORElem
Compute weak normal form of `f` with respect to `G` and return the remainder `r` and quotients `q`.
"""
function wnormalform(f::T, G::Vector{T}; order::Symbol=:grevlex) where T<: DORElem
    r = f
    n = nvars(parent(f))
    q = [zero(parent(f)) for _ in G]

    reducible = true
    while reducible 
        reducible = false
        for (i, g) in enumerate(G)
            r == zero(parent(f)) && break

            lt_r = leading_term(r, order=order)
            lt_g = leading_term(g, order=order)
            lc_r = coefficients(lt_r)[1]
            lc_g = coefficients(lt_g)[1]
            a = exponent_vectors(lt_r)[1] - exponent_vectors(lt_g)[1]

            if minimum(a) >= 0
                reducible = true
                q_mon = one(parent(f))
                for j = 1:n
                    q_mon *= dgens(parent(f))[j] ^ a[j]
                end
                q[i] = DORElem(parent(f), unwrap(q[i]) + lc_r // lc_g  * unwrap(q_mon))
                r = r - DORElem(parent(f), lc_r // lc_g  * unwrap(q_mon)) * g
                break
            end
        end
    end
    return r, q 
    
end

"""
    pfaffian_system(G::Vector{T}, S::Vector{T}) where T <: DORElem

# Examples

```jldoctest
julia> R, (x, y), (dx, dy) = diff_op_ring(["x", "y"])
(2-dimensional ring of differential opeartors in [x,y], PfaffianSystems.DORElem{AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.RationalFunctionFieldElem{Rational{BigInt}, AbstractAlgebra.Generic.MPoly{Rational{BigInt}}}}}[x, y], PfaffianSystems.DORElem{AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.RationalFunctionFieldElem{Rational{BigInt}, AbstractAlgebra.Generic.MPoly{Rational{BigInt}}}}}[dx, dy])
julia> pfaffian_system([dx^2 + 1, dy^2 + 1], [one(dx), dx, dy])
2-element Vector{Matrix{AbstractAlgebra.Generic.RationalFunctionFieldElem{Rational{BigInt}, AbstractAlgebra.Generic.MPoly{Rational{BigInt}}}}}:
 [0 1 0; -1 0 0; 0 0 0]
 [0 0 1; 0 0 0; -1 0 0]
```
"""
function pfaffian_system(G::Vector{T}, S::Vector{T}) where T <: DORElem

    R = base_ring(parent(S[1]))
    dops = dgens(parent(S[1]))
    n = size(dops)[1]
    d = size(S)[1]

    # p = Vector{Matrix{elem_type(R)}}(Matrix{elem_type(R)}(undef, d, d), n)
    p = [fill(zero(R), d, d) for _ in 1:n]
    for i in eachindex(p)
        for j = 1:d
            nf_ds1 = normalform(dops[i] * S[j], G, order=:grevlex)[1]
            nf_ds1_coef = coefficients(nf_ds1)
            nf_ds1_mono = monomials(nf_ds1)

            for k = 1:d
                idx = findfirst(x->x==S[k], nf_ds1_mono)
                if idx === nothing
                    p[i][j, k] = zero(R)
                else
                    p[i][j, k] = nf_ds1_coef[idx]
                end
            end
        end
    end

    return p
end
