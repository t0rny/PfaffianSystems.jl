"""
	addVars(name::AbstractString[, n::Integer, v2d::Bijection])

Generate a Symbolics variable `name` and its corresponding differential operator d`name`. When `n` is given, generate vectors [`name`1...`name`n] and [d`name`1...d`name`n]. 
The correspondence between `name` and d`name` is added to the bijection `v2d`. 

# Examples

```jldoctest
julia> x, dx, v2d = genVars("x")
(x, dx, Bijection{Symbolics.Num,Symbolics.Num} (with 1 pairs))

julia> y, dy, v2d = addVars("y", 2, v2d)
(Symbolics.Num[y1, y2], Symbolics.Num[dy1, dy2], Bijection{Symbolics.Num,Symbolics.Num} (with 3 pairs))
```
"""
function addVars(name::AbstractString, v2d::Bijection{Num, Num})
	var_ex = Symbol(name)
	diffop_ex = Symbol("d", var_ex)
	var, diffop = @variables $var_ex::Rational, $diffop_ex::Rational
	v2d[var] = diffop
	return var, diffop, v2d
end

function addVars(name::AbstractString, n::Integer, v2d::Bijection{Num, Num})
	vars = Vector{Num}(undef, n)
	diffops = Vector{Num}(undef, n)

	for i = 1:n
		vars[i], diffops[i], v2d = addVars(name*string(i), v2d)
	end

	return vars, diffops, v2d
end

"""
	genVars(name::AbstractString[, n::Integer])

Generate a Symbolics variable `name` and its corresponding differential operator d`name`. When `n` is given, generate vectors [`name`1...`name`n] and [d`name`1...d`name`n]. 
The correspondence between `name` and d`name` is returned as a new bijection `v2d`, which will be provided to `addVars`. 

# Examples

```jldoctest
julia> x, dx, v2d = genVars("x")
(x, dx, Bijection{Symbolics.Num,Symbolics.Num} (with 1 pairs))

julia> y, dy, v2d = addVars("y", 2, v2d)
(Symbolics.Num[y1, y2], Symbolics.Num[dy1, dy2], Bijection{Symbolics.Num,Symbolics.Num} (with 3 pairs))
```
"""
genVars(name::AbstractString) = addVars(name, Bijection{Num, Num}())
genVars(name::AbstractString, n::Integer) = addVars(name, n, Bijection{Num, Num}())

function asir_derivative(sym::Num, var::Num)
	vars_list = get_variables(sym) .|> Num
	asir_cmd = "diff($sym, $var);"
	# """
	# diff($sym, $var);
	# """
	asir_res = asir_cmd |> runAsir |> parseAsir
	return evalAsir(asir_res[1], vars_list)
end

"""
	apply_dmon(DOmon::AbstractTerm, F::Num, p2s::Bijection, v2d::Bijection)

Apply a term of differential operator `DOterm` to an expression `F`. 
The bijection `p2s` relates a MulltivariatePolunomials expression to a SymbolicUtils one, and `d2v` does a differential operator to its corresponding variable. 
"""
function apply_doterm(DOterm::AbstractTerm, F::Num, p2s::Bijection, v2d::Bijection{Num, Num})
	d2v = inv(v2d)
	coef = coefficient(DOterm)
	mon = monomial(DOterm)
	diffops = variables(mon) .|> (s->p2s[s])
	exps = exponents(mon)

	retF = F
	for (diffop, e) in zip(diffops, exps)
		if haskey(d2v, diffop)
			for k = 1:e
				retF = asir_derivative(retF, d2v[diffop])
			end
		else
			coef = diffop^e*coef
		end
	end
	return coef*retF
end
function apply_doterm(DOterm::PolyForm, F::Num, v2d::Bijection{Num, Num})
	@assert length(terms(DOterm.p)) == 1 "Error: Differnetial operator has more than 1 term"
	apply_doterm(terms(DOterm.p)[1], F, DOterm.pvar2sym, v2d)
end
apply_doterm(DOterm::Num, F::Num, v2d::Bijection{Num, Num}) = apply_doterm(DOterm |> value |> PolyForm, F, v2d)


function apply_do(DiffOp::BasicSymbolic, F::Num, v2d::Bijection{Num, Num})
	dp = PolyForm(DiffOp)
	p2s = dp.pvar2sym
	retF::Num = 0
	for DOterm in terms(dp.p)
		retF = retF + apply_doterm(DOterm, F, p2s, v2d)
	end
	return retF
end
function apply_do(DiffOp::Union{Integer, AbstractFloat}, F::Num, v2d::Bijection{Num, Num})
	return DiffOp*F
end
apply_do(DiffOp::Num, F::Num, v2d::Bijection) = apply_do(DiffOp |> value, F, v2d)