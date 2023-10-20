
############################################################
# 
# Ideal of Weyl algebra
# 
############################################################

mutable struct DIdeal{T <: AbstractDiffOp} <: Ideal{T}
    base_ring::AbstractDORing
    gens::Vector{T}

    function DIdeal{T}(R::AbstractDORing, gens::Vector) where T <: AbstractDiffOp
        if !all([R == parent(g) for g in gens])
            # error("All generators must be elements of the same Weyl algebra")
            throw(ArgumentError("All generators must be included in", R))
        end
        new{T}(R, gens)
    end 
end

# function DIdeal(R::AbstractDORing, gens::Vector{T}) where T <: AbstractDiffOp
#     DIdeal{eltype(gens)}(R, gens)
# end
function DIdeal(R::AbstractDORing, gens::Vector)
		DIdeal{elem_type(R)}(R, R.(gens))
end
DIdeal(gens::Vector{T}) where T <: AbstractDiffOp = DIdeal(parent(gens[1]), gens)


base_ring(I::DIdeal) = I.base_ring
gens(I::DIdeal) = I.gens

function Base.show(io::IO, I::DIdeal)
    print(io, "Ideal of ", I.base_ring, " generated by [", join(string.(I.gens), ", "), "]")
end

function (D::AbstractDORing)(I::DIdeal)
    DIdeal(D, I |> gens .|> (s->coerce(s, D)))
end

"""
	evaluate(I::DIdeal{T}, vrs::Vector{T}, vals::Vector{T}) where T <: WAlgElem

Return the ideal that is obtained by evaluating `I` at `vrs` equal to `vals`.

#Examples
```jldoctest
julia> D, (x, y), (dx, dy) = weyl_algebra(["x", "y"])
(2-d Weyl algebra in [x,y], PfaffianSystems.WAlgElem{AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.MPoly{Rational{BigInt}}}}[x, y], PfaffianSystems.WAlgElem{AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.MPoly{Rational{BigInt}}}}[dx, dy])
julia> I = DIdeal(D, [dx + x])
Ideal of 2-d Weyl algebra in [x,y] generated by [dx + x]
julia> evaluate(I, [x, dx], [y, dy])
Ideal of 2-d Weyl algebra in [x,y] generated by [dy + y]
```
"""
evaluate(I::DIdeal{T}, vrs::Vector{T}, vals::Vector{T}) where T <: WAlgElem = DIdeal(base_ring(I), evaluate.(gens(I), (vrs,), (vals,)))


"""
	intersection_DIdeals(Is::Ideal...)

Return the intersection of the D-ideals `Is[1]`, `Is[2]`, ...

# Examples

```jldoctest
julia> D, x, dx = weyl_algebra("x")
(1-d Weyl algebra in [x], x, dx)
julia> I1 = DIdeal(D, [dx^2 + 1])
Ideal of 1-d Weyl algebra in [x] generated by [dx^2 + 1]
julia> I2 = DIdeal(D, [dx + x])
Ideal of 1-d Weyl algebra in [x] generated by [dx + x]
julia> intersection_DIdeals(I1, I2)
Ideal of 1-d Weyl algebra in [x] generated by [x*dx^3 + (x^2 - 2)*dx^2 + x*dx + x^2 - 2, dx^5 + x*dx^4 + 6*dx^3 + 2*x*dx^2 + 5*dx + x]
```
"""
function intersection_DIdeals(Is::DIdeal{T}...) where T <: WAlgElem
    D = base_ring(Is[1])
    m = length(Is)
    all([D == base_ring(I) for I in Is]) || throw(DomainError("Base rings must be the same"))

    Dt, vall, dvall = weyl_algebra(D, ["t"*string(i) for i in 1:m]; new_vars_pos=:prepend)
		t = vall[1:m]
		v = vall[m+1:end]
    # dv = dvall[1:end-m]
    # t = vall[end-m+1:end]
		dt = dvall[1:m]

    genJ = [reduce(+, t)-1]
    for (i, I) in enumerate(Is)
        append!(genJ, [t[i]*g for g in Dt.(gens(I))])
    end

	asir_cmd = 
	"""
	load("nk_restriction.rr")\$
	V = [$(vec2str([vall; dvall]))]\$
	M = nk_restriction.make_elim_order_mat($m, $(length(v)))\$
	J = [$(vec2str(genJ))]\$
	GB = nd_weyl_gr(J, V, 0, M);
	"""
	asir_res = asir_cmd |> runAsir |> parseAsir
	asir_res = filter!((s)->(startswith(s, "[")), asir_res)

	(length(asir_res) != 1) && throw(DomainError("Invalid result from Asir", asir_res))

	elimOrdGens = Dt.(evalAsir(asir_res[1], [vall; dvall]))

	notHasTmpVar(dop) = (vars(dop), [t; dt]) |> (s-> intersect(s...)) |> isempty

	return D(DIdeal(filter(notHasTmpVar, elimOrdGens)))
end

"""
	restriction_DIdeal(I::DIdeal{T}, rest_vars::OrderedSet{T}) where T <: WAlgElem
	restriction_DIdeal(I::DIdeal{T}, rest_vars::Vector{T}) where T <: WAlgElem

Return the restriction ideal of `I` with respect to the subspace where every variable of `rest_vars` equals to zero.

# Examples
```jldoctest
julia> D, (x, y), (dx, dy) = weyl_algebra(["x", "y"])
(2-d Weyl algebra in [x,y], PfaffianSystems.WAlgElem{AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.MPoly{Rational{BigInt}}}}[x, y], PfaffianSystems.WAlgElem{AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.MPoly{Rational{BigInt}}}}[dx, dy])
julia> I = DIdeal(D, [dx^2 + 2*dx*dy + dy^2 + 1, -dy + x - y])
Ideal of 2-d Weyl algebra in [x,y] generated by [dx^2 + 2*dx*dy + dy^2 + 1, -dy + x - y]
julia> restriction_DIdeal(I, [y])
Ideal of 2-d Weyl algebra in [x,y] generated by [-dx^2 - 2*x*dx - x^2 - 2]
```
"""
function restriction_DIdeal(I::DIdeal{T}, rest_vars::OrderedSet{T}) where T <: WAlgElem
	D = base_ring(I)
	vs = gens(D)
	dos = dgens(D)
	n = length(vs)
	!issubset(rest_vars, vs) && throw(DomainError("'unknown variables: $(setdiff(rest_vars, vs)) in the second argument'"))

	res_indcs = [findfirst((s->s == v), vs) for v in rest_vars]
	res_vs = vs[res_indcs]
	res_dos = dos[res_indcs]

	rem_indcs = setdiff(1:n, res_indcs)
	rem_vs = vs[rem_indcs]
	rem_dos = dos[rem_indcs]
	m = length(res_vs)

	asir_cmd = 
	"""
	load("nk_restriction.rr")\$
	V = [$(vec2str(res_vs, rem_vs))]\$
	DV = [$(vec2str(res_dos, rem_dos))]\$
	W = [$(vec2str((1:n .< m+1) .|> Int |> collect))]\$
	J = [$(vec2str(gens(I)))]\$
	GB = nk_restriction.restriction_ideal(J, V, DV, W);
	GB;
	"""

	asir_res = asir_cmd |> runAsir |> parseAsir
	asir_res = filter!((s)->(startswith(s, "[")), asir_res)

	(length(asir_res) != 1) && throw(DomainError("Invalid result from Asir: $asir_res"))

	vars_list = cat(rem_vs, rem_dos; dims=1)
	restGens = D.(evalAsir(asir_res[1], vars_list))

	return DIdeal(base_ring(I), restGens)
end
restriction_DIdeal(I::DIdeal{T}, rest_vars::Vector{T}) where T <: WAlgElem = restriction_DIdeal(I, OrderedSet(rest_vars))

"""
	integration_DIdeal(I::DIdeal{T}, integ_vars::OrderedSet{T}) where T <: WAlgElem
	integration_DIdeal(I::DIdeal{T}, integ_vars::Vector{T}) where T <: WAlgElem

Return the integration ideal of `I` with respect to `integ_vars`.

# Examples
```jldoctest
julia> D, (x, y), (dx, dy) = weyl_algebra(["x", "y"])
(2-d Weyl algebra in [x,y], PfaffianSystems.WAlgElem{AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.MPoly{Rational{BigInt}}}}[x, y], PfaffianSystems.WAlgElem{AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.MPoly{Rational{BigInt}}}}[dx, dy])
julia> I = DIdeal(D, [dx+x, dy+y])
Ideal of 2-d Weyl algebra in [x,y] generated by [dx + x, dy + y]
julia> integration_DIdeal(I, [x])
Ideal of 2-d Weyl algebra in [x,y] generated by [dy + y]
```
"""
function integration_DIdeal(I::DIdeal{T}, integ_vars::OrderedSet{T}) where T <: WAlgElem
	D = base_ring(I)
	vs = gens(D)
	dos = dgens(D)
	n = length(vs)
	!issubset(integ_vars, vs) && throw(DomainError("'unknown variables: $(setdiff(integ_vars, vs)) in the second argument'"))

	int_indcs = [findfirst((s->s == v), vs) for v in integ_vars]
	int_vs = vs[int_indcs]
	int_dos = dos[int_indcs]

	rem_indcs = setdiff(1:n, int_indcs)
	rem_vs = vs[rem_indcs]
	rem_dos = dos[rem_indcs]
	m = length(int_vs)

	asir_cmd = 
	"""
	load("nk_restriction.rr")\$
	V = [$(vec2str(int_vs, rem_vs))]\$
	DV = [$(vec2str(int_dos, rem_dos))]\$
	W = [$(vec2str((1:n .< m+1) .|> Int |> collect))]\$
	J = [$(vec2str(gens(I)))]\$
	GB = nk_restriction.integration_ideal(J, V, DV, W);
	GB;
	"""

	asir_res = asir_cmd |> runAsir |> parseAsir
	asir_res = filter!((s)->(startswith(s, "[")), asir_res)

	(length(asir_res) != 1) && throw(DomainError("Invalid result from Asir", asir_res))

	vars_list = cat(rem_vs, rem_dos; dims=1)
	integGens = D.(evalAsir(asir_res[1], vars_list))

	return DIdeal(base_ring(I), integGens)
end
integration_DIdeal(I::DIdeal{T}, integ_vars::Vector{T}) where T <: WAlgElem = integration_DIdeal(I, OrderedSet(integ_vars))

"""
	multiplication_DIdeal(I::DIdeal{T}, J::DIdeal{T}) where T <: WAlgElem

Return an ideal that annihilates the product `fg` where `f` and `g` are annihilated by `I` and `J`, respectively. 

# Exmaples
```jldoctest
julia> D, x, dx = weyl_algebra("x")
(1-d Weyl algebra in [x], x, dx)
julia> I1 = DIdeal(D, [dx^2 + 1])
Ideal of 1-d Weyl algebra in [x] generated by [dx^2 + 1]
julia> I2 = DIdeal(D, [dx + x])
Ideal of 1-d Weyl algebra in [x] generated by [dx + x]
julia> multiplication_DIdeal(I1, I2)
Ideal of 1-d Weyl algebra in [x] generated by [-dx^2 - 2*x*dx - x^2 - 2]
```
"""
function multiplication_DIdeal(I::DIdeal{T}, J::DIdeal{T}) where T <: WAlgElem
	D = base_ring(I)
	D2, vall, dvall = weyl_algebra(D, ["t"*string(i) for i in 1:nvars(D)])

	vs = vall[1:nvars(D)]
	dos = dvall[1:nvars(D)]
	vs_dummy = vall[nvars(D)+1:end]
	dos_dummy = dvall[nvars(D)+1:end]

	Jgens = [evaluate(g, [vs; dos], [vs_dummy; dos_dummy]) for g in gens(D2(J))]

	Igens = [evaluate(g, dos, dos .+ dos_dummy) for g in gens(D2(I))]
	Jgens = [evaluate(g, [vs_dummy; dos_dummy], [vs .- vs_dummy; -dos_dummy]) for g in Jgens]

	return restriction_DIdeal(DIdeal(D2, [Igens; Jgens]), vs_dummy) |> D
end

function mult_test()
D, (x, s1, s2, s3), (dx, ds1, ds2, ds3) = weyl_algebra(["x", "s1", "s2", "s3"])

IHdiff1 = DIdeal(D, [dx^2 - 2*x*dx + 2*(-ds1*s1-1), (-ds1*s1+1)*s1^2 - 2*x*(-ds1*s1+2)*s1 + 2*(-ds1*s1+2)*(-ds1*s1+1), ds2, ds3])

Igauss = DIdeal(D, [dx + 2*x, ds1, ds2, ds3])

multiplication_DIdeal(IHdiff1, Igauss)

end

"""
	stdmon!(I::DIdeal[, ordered_vars::Vector{Num}])

Compute standard monomials if `I` is zero-dimensional as an ideal of R-ring. The order of varaibles in the computation of Grobner bases can be provided as `ordered_vars`, which can be omitted as long as the computational time is not crutial. 
"""
function stdmon(I::DIdeal{T}, ordered_vars::OrderedSet{T}) where T <: DORElem
	asir_cmd = 
	"""
	load("yang.rr")\$
	V = [$(ordered_vars |> collect |> vec2str)]\$
	yang.define_ring(["partial", V])\$
	Gb = yang.gr([$(vec2str(gens(I)))])\$
	In = yang.in(Gb)\$
	yang.stdmon(Gb);
	"""

	asir_res = runAsir(asir_cmd) |> parseAsir
	asir_res = filter!((s->startswith(s, "[")), asir_res)

	(length(asir_res) != 1) && throw(DomainError("Invalid result from Asir", asir_res))

	D = base_ring(I)
	vars_list = [gens(D); dgens(D)]
	return D.(evalAsir(asir_res[1], vars_list)) |> reverse
end
stdmon(I::DIdeal) = stdmon(I, OrderedSet(gens(base_ring(I))))
stdmon(I::DIdeal, vars::Vector) = stdmon(I, vars |> OrderedSet)

"""
    pfaffian_system2(G::Vector{T}, S::Vector{T}) where T <: DORElem

# Examples

```jldoctest
julia> R, (x, y), (dx, dy) = diff_op_ring(["x", "y"])
(2-dimensional ring of differential opeartors in [x,y], PfaffianSystems.DORElem{AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.RationalFunctionFieldElem{Rational{BigInt}, AbstractAlgebra.Generic.MPoly{Rational{BigInt}}}}}[x, y], PfaffianSystems.DORElem{AbstractAlgebra.Generic.MPoly{AbstractAlgebra.Generic.RationalFunctionFieldElem{Rational{BigInt}, AbstractAlgebra.Generic.MPoly{Rational{BigInt}}}}}[dx, dy])
julia> pfaffian_system2([dx^2 + 1, dy^2 + 1], [one(dx), dx, dy])
2-element Vector{Matrix{AbstractAlgebra.Generic.RationalFunctionFieldElem{Rational{BigInt}, AbstractAlgebra.Generic.MPoly{Rational{BigInt}}}}}:
 [0 1 0; -1 0 0; 0 0 0]
 [0 0 1; 0 0 0; -1 0 0]
```
"""
function pfaffian_system2(I::DIdeal{T}, sm::OrderedSet{T}) where T <: DORElem
	v = gens(base_ring(I))
	# v = collect(ordered_vars)
	n = length(v)
	# sm = stdmon(I, ordered_vars)
	d = length(sm)
	asir_cmd = 
	"""
	load("yang.rr")\$
	V = [$(vec2str(v))]\$
	yang.define_ring(["partial", V])\$
	Gb = yang.gr([$(vec2str(gens(I)))])\$
	Pf = yang.pf([$(vec2str(sm))], Gb)\$
	for (I=0; I < $n; I++){
		print(Pf[I])\$
	}
	"""

	@show asir_cmd

	asir_res = runAsir(asir_cmd) |> parseAsir
	asir_res = filter!((s->startswith(s, "[")), asir_res)

	size(asir_res, 1) != length(v)*length(sm) && throw(DomainError("Invalid result from Asir", asir_res))

	R = base_ring(I) |> coef_ring
	p = [fill(zero(R), d, d) for _ = 1:n]

	for i in eachindex(p)
		for j = 1:d
			p[i][j, :] = R.(evalAsir(asir_res[(i-1)*d + j], gens(R)))
		end
	end

	# (length(asir_res) != 1) && throw(DomainError("Invalid result from Asir", asir_res))
	return p
end
pfaffian_system2(I::DIdeal{T}, vars::Vector{T}) where T <: DORElem = pfaffian_system2(I, OrderedSet(vars))
pfaffian_system2(I::DIdeal) = pfaffian_system2(I, OrderedSet(gens(base_ring(I))))

# function stdmon!(I::DIdeal, ordered_vars::OrderedSet{Num})
# 	@assert isequal(I.v2d.domain, ordered_vars) "Error: invalid vector for second argument"
# 	asir_cmd = 
# 	"""
# 	load("yang.rr")\$
# 	V = [$(vec2str(ordered_vars))]\$
# 	yang.define_ring(["partial", V])\$
# 	Gb = yang.gr([$(vec2str(I.gens))])\$
# 	In = yang.in(Gb)\$
# 	yang.stdmon(Gb);
# 	"""

# 	asir_res = runAsir(asir_cmd) |> parseAsir
# 	asir_res = filter!((s->startswith(s, "[")), asir_res)

# 	# (length(asir_res) != 1) && return nothing
# 	if length(asir_res) != 1
# 		I.flags["isZeroDim"] = false
# 		return Vector{Num}()
# 	end

# 	I.flags["isZeroDim"] = true
# 	vars_list = cat(collect(I.v2d.domain), collect(I.v2d.range); dims=1)
# 	return evalAsir(asir_res[1], vars_list)
# end
# stdmon!(I::DIdeal) = stdmon!(I, I.v2d.domain |> OrderedSet{Num})
# stdmon!(I::DIdeal, vars::Vector{Num}) = stdmon!(I, vars |> OrderedSet{Num})


# """
# 	DIdeal <: AbstractIdeal

# D-Ideal type
# """
# struct DIdeal <: AbstractIdeal
# 	gens::Vector{Num} 
# 	v2d::Bijection{Num, Num}
# 	flags::Dict{String, Union{Nothing, Bool}}

# 	function DIdeal(gens::Vector{Num}, v2d::Bijection{Num, Num})
# 		gen_vars = get_variables.(gens) .|> Set |> (s->union(s...))
# 		dict_vars = union(v2d.domain, v2d.range)

# 		@assert issubset(gen_vars, dict_vars) "Error: generators include unknown variables"
# 		_v2d = empty(v2d)
# 		for (k, v) in v2d
# 			if k in gen_vars || v in gen_vars
# 				_v2d[k] = v
# 			end
# 		end
# 		new(copy(gens), _v2d, Dict{String, Union{Nothing, Bool}}("isZeroDim"=>nothing))
# 	end
# 	DIdeal(gens::AbstractVector, v2d::Bijection{Num, Num}) = DIdeal(Num.(gens), v2d)
# end

# function makeTestVarsAndIdeal()
# 	x, dx, v2d = genVars("x", 3)
# 	# return x, dx, v2d, DIdeal([dx[1]^2 + 1, x[2]*dx[2] - 2, x[3]*dx[3] - 1], v2d)
# 	return x, dx, v2d, DIdeal([dx[1] + 2*x[1], dx[2]^2 + 1, x[3]*dx[3] - 1], v2d)
# end

# function apply_ideal(I::DIdeal, F::Num; use_asir=false)
# 	map(I.gens) do DO
# 		apply_do(DO, F, I.v2d; use_asir=use_asir)
# 	end
# end




# """
# 	isZeroDimensional(I::DIdeal)

# Check if the ideal in R generated by DIdeal `I` is zero-dimensional. 
# """
# function isZeroDimensional(I::DIdeal; std_mons=Vector())
# 	if !isnothing(I.flags["isZeroDim"])
# 		return I.flags["isZeroDim"]
# 	end
# 	stdmon!(I) |> (s->append!(std_mons, s))
# 	return I.flags["isZeroDim"]
# end


# """
# 	eliminationIdeal(I::DIdeal, elim_vars::OrderedSet{Num})
# 	eliminationIdeal(I::DIdeal, elim_vars::Vector{Num})

# Return elimination ideal of `I` that does not include any variable of `elim_vars`. 
# """
# function eliminationIdeal(I::DIdeal, elim_vars::OrderedSet{Num})
# 	@assert issubset(elim_vars, I.v2d.domain) "Error: unknown variables in the second argument"

# 	elim_dos = map(s->I.v2d[s], elim_vars |> collect)
# 	rem_vars = setdiff(I.v2d.domain, elim_vars) |> collect
# 	rem_dos = map(s->I.v2d[s], rem_vars |> collect)
# 	elim_vars = collect(elim_vars)
# 	m = length(elim_vars)
# 	n = length(I.v2d.domain)

# 	asir_cmd = 
# 	"""
# 	load("nk_restriction.rr")\$
# 	V = [$(vec2str(elim_vars, rem_vars, elim_dos, rem_dos))]\$
# 	M = nk_restriction.make_elim_order_mat($(2*m), $(2*(n - m)))\$
# 	J = [$(vec2str(I.gens))]\$
# 	GB = nd_weyl_gr(J, V, 0, M);
# 	"""

# 	asir_res = asir_cmd |> runAsir |> parseAsir
# 	asir_res = filter!((s)->(startswith(s, "[")), asir_res)

# 	(length(asir_res) != 1) && return nothing

# 	vars_list = cat(rem_vars, elim_vars, rem_dos, elim_dos; dims=1)
# 	elimOrdGens = evalAsir(asir_res[1], vars_list)

# 	notHasElimVar(diffOp) = (get_variables(diffOp), [elim_vars; elim_dos]) .|> Set |> (s->intersect(s...)) |> isempty

# 	return DIdeal(filter(notHasElimVar, elimOrdGens), Dict((v, d) for (v,d) in zip(rem_vars, rem_dos)) |> Bijection)
# end
# eliminationIdeal(I::DIdeal, elim_vars::Vector{Num}) = eliminationIdeal(I, OrderedSet(elim_vars))



# function restrictionIdeal(I::DIdeal, rest_vars::OrderedSet{Num})
# 	@assert issubset(rest_vars, I.v2d.domain) "Error: unknown variables in the second argument"

# 	rest_vs = rest_vars |> collect |> sort
# 	rem_vs = setdiff(I.v2d.domain, rest_vars) |> collect |> sort

# 	rest_dos = map(s->I.v2d[s], rest_vs)
# 	rem_dos = map(s->I.v2d[s], rem_vs)

# 	m = length(rest_vs)
# 	n = m + length(rem_vs)

# 	asir_cmd = 
# 	"""
# 	load("nk_restriction.rr")\$
# 	V = [$(vec2str(rest_vs, rem_vs))]\$
# 	DV = [$(vec2str(rest_dos, rem_dos))]\$
# 	M = [$(vec2str((1:n .< m+1) .|> Int |> collect))]\$
# 	J = [$(vec2str(I.gens))]\$
# 	GB = nk_restriction.restriction_ideal(J, V, DV, M);
# 	GB;
# 	"""

# 	asir_res = asir_cmd |> runAsir |> parseAsir
# 	asir_res = filter!((s)->(startswith(s, "[")), asir_res)

# 	(length(asir_res) != 1) && return nothing

# 	vars_list = cat(rem_vs, rem_dos; dims=1)
# 	restGens = evalAsir(asir_res[1], vars_list)

# 	return DIdeal(restGens, ((v, d) for (v,d) in zip(rem_vs, rem_dos)) |> Dict |> Bijection)
# end
# restrictionIdeal(I::DIdeal, rest_vars::Vector{Num}) = restrictionIdeal(I, OrderedSet(rest_vars))
