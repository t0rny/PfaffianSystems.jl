
# PolyDiffOp(a::Vector{T}, x::MonomialVector{false}) = Polynomial{false, Rational{Integer}}(a, x)
const DGen = PolyVar{false}
abstract type DiffOp end

# struct Dgen <: DiffOp
# 	p::PolyVar{false}
# 	v2d::Bijection{PolyVar{false}, PolyVar{false}}
# 	p2s::Bijection{PolyVar{false}, Num}
# end

struct PolyDiffOp <: DiffOp
	p::Polynomial{false, Rational}
	v2d::Bijection{DGen, DGen}
	p2s::Bijection{DGen, Num}

	PolyDiffOp(p::Polynomial) = new(p, get_var2diff(), get_pv2sym())
	# PolyDiffOp(g::DGen) = new(convert(Polynomial{false, Rational}, g), get_var2diff(), get_pv2sym())
	PolyDiffOp(g::DGen, v2d, p2s) = new(convert(Polynomial{false, Rational}, g), v2d, p2s)
end

Base.print(io::IO, dop::DiffOp) = print(io, dop.p)
Base.show(io::IO, dop::DiffOp) = show(io, dop.p)
# const PolyDiffOp = Polynomial{false, Rational}
# PolyDiffOp(p::Polynomial) = convert(PolyDiffOp, p)

# Refer to construction of PolyForm in [SymbolicUtils](https://github.com/JuliaSymbolics/SymbolicUtils.jl)
const VAR2DIFF = Ref(WeakRef(nothing))
const PV2SYM = Ref(WeakRef(nothing))

clear_dicts() = begin
	VAR2DIFF[] = WeakRef(nothing)
	PV2SYM[] = WeakRef(nothing)
	nothing
end

function get_var2diff()
	v = VAR2DIFF[].value
	if v === nothing
		d = Bijection{DGen, DGen}()
		VAR2DIFF[] = WeakRef(d)
		return d
	else
		return v
	end
end

function get_pv2sym()
	v = PV2SYM[].value
	if v === nothing
		d = Bijection{DGen, Num}()
		PV2SYM[] = WeakRef(d)
		return d
	else
		return v
	end
end

# function check_dicts_validity(p, q)
# 	@assert p.v2d === q.v2d "v2d mismatch"
# 	@assert p.p2s === q.p2s "v2d mismatch"

# 	return p.v2d, p.p2s
# end

Base.:(+)(p::PolyDiffOp, q::PolyDiffOp) = PolyDiffOp(p.p + q.p)

Base.:(*)(p::PolyDiffOp, q::PolyDiffOp) = PolyDiffOp(p.p * q.p)

Base.:(^)(p::PolyDiffOp, n::Integer) = PolyDiffOp(p.p^n)

# TOOD: differential operators with coefficients in rational functions should be implemented
# const RatDiffOp = Polynomial{false, RationalPoly{Polynomial{true, Rational{T}}, Polynomial{true, Rational{T}}}} where T
# function RatDiffOp(p::Polynomial{C, T}) where {C, T}
# end

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

function genVars2(name::AbstractString, v2d=get_var2diff(), p2s=get_pv2sym())
	# v2d = get_var2diff()
	# p2s = get_pv2sym()

	@assert !startswith(name, 'd') "variable name \"$name\" must start from letters except for \"d\""
	# @assert !in(name, DP.name.(v2d.domain)) "already exists variable named \"$name\""

	var_ex = Symbol(name)
	dop_ex = Symbol('d', var_ex)
	symvar, symdop = @variables $var_ex, $dop_ex

	s2p = active_inv(p2s)

	if haskey(s2p, symvar) && haskey(s2p, symdop) && v2d[s2p[symvar]] == s2p[symdop]
		@show objectid(get_var2diff())
		@show objectid(v2d)
		@show objectid(get_pv2sym())
		@show objectid(p2s)
		@warn "variable already exist"
		return PolyDiffOp(s2p[symvar], v2d, p2s), PolyDiffOp(s2p[symdop], v2d, p2s)
	else

		pvar = DGen(name)
		pdop = DGen('d'*name)

		s2p[symvar] = pvar
		s2p[symdop] = pdop
		v2d[pvar] = pdop
		@show objectid(get_var2diff())
		@show objectid(v2d)
		@show objectid(get_pv2sym())
		@show objectid(p2s)
		return PolyDiffOp(pvar, v2d, p2s), PolyDiffOp(pdop, v2d, p2s)
	end
end
# genVars2(name::AbstractString) = genVars2(name, get_var2diff(), get_pv2sym())

function genVars2(name::AbstractString, n::Integer)
	pvars = Vector{DiffOp}(undef, n)
	pdiffops = Vector{DiffOp}(undef, n)

	for i in eachindex(pvars)
		pvars[i], pdiffops[i] = genVars2(name*string(i))
	end

	return pvars, pdiffops
end

# genVars2(name::AbstractString) = addVars2(name)
# genVars2(name::AbstractString) = addVars2(name, Bijection{Num, Num}())


"""
	apply_dmon(DOmon::AbstractTerm, F::Num, p2s::Bijection, v2d::Bijection)

Apply a term of differential operator `DOterm` to an expression `F`. 
The bijection `p2s` relates a MulltivariatePolunomials expression to a SymbolicUtils one, and `d2v` does a differential operator to its corresponding variable. 
"""
function apply_doterm(DOterm::AbstractTerm, F::Num, p2s::Bijection, v2d::Bijection{Num, Num}; use_asir=false)
	d2v = inv(v2d)
	coef = coefficient(DOterm)
	mon = monomial(DOterm)
	diffops = variables(mon) .|> (s->p2s[s])
	exps = exponents(mon)

	retF = F
	for (diffop, e) in zip(diffops, exps)
		if haskey(d2v, diffop)
			for k = 1:e
				retF = use_asir ? asir_derivative(retF, d2v[diffop]) : derivative(retF, d2v[diffop])
			end
		else
			coef = diffop^e*coef
		end
	end
	return coef*retF
end
function apply_doterm(DOterm::PolyForm, F::Num, v2d::Bijection{Num, Num}; use_asir=false)
	@assert length(terms(DOterm.p)) == 1 "Error: Differnetial operator has more than 1 term"
	apply_doterm(terms(DOterm.p)[1], F, DOterm.pvar2sym, v2d; use_asir=use_asir)
end
apply_doterm(DOterm::Num, F::Num, v2d::Bijection{Num, Num}; use_asir=false) = apply_doterm(DOterm |> value |> PolyForm, F, v2d; use_asir=use_asir)


function apply_do(DiffOp::BasicSymbolic, F::Num, v2d::Bijection{Num, Num}; use_asir=false)
	dp = PolyForm(DiffOp)
	p2s = dp.pvar2sym
	retF::Num = 0
	for DOterm in terms(dp.p)
		# retF = retF + apply_doterm(DOterm, F, p2s, v2d; use_asir=use_asir)
		retF += apply_doterm(DOterm, F, p2s, v2d; use_asir=use_asir)
	end
	return retF
end
function apply_do(DiffOp::Union{Integer, AbstractFloat}, F::Num, v2d::Bijection{Num, Num}; use_asir=false)
	return DiffOp*F
end
apply_do(DiffOp::Num, F::Num, v2d::Bijection{Num, Num}; use_asir=false) = apply_do(DiffOp |> value, F, v2d; use_asir=use_asir)

# function apply_do2(DiffOp::PolyDiffOp, F::Num, v2d::{})
# end

isDiff(s::DGen) = startswith(DP.name(s), 'd')
function canonic(diffop::PolyDiffOp)
	v2d = diff
	for do_term in terms(diffop)
		mon = monomial(do_term)
		vars = variables(mon)
		exps = exponents(mon)

		for i in reverse(eachindex(exps))
		end
	end
end

function dmul(dol::Num, dor::BasicSymbolic, v2d::Bijection{Num, Num}; use_asir=false)
	dor_pf = PolyForm(dor)
	p2s = dor_pf.pvar2sym
	retDO = dor*dol
	for dor_term in terms(dor_pf.p)
		# retDO += Num(dor_term)*dol
		coef = coefficient(dor_term)
		mon = monomial(dor_term)
		syms = variables(mon) .|> (s->p2s[s])
		exps = exponents(mon)
		varIdx = map(syms) do s s in v2d.domain end
		var_term = reduce(*, [Num(s^e) for (s, e) in zip(syms[varIdx], exps[varIdx])])

		retDO += coef*apply_do(dol, Num(var_term), v2d; use_asir=use_asir)*reduce(*, [Num(s^e) for (s, e) in zip(syms[.!varIdx], exps[.!varIdx])])
	end
	return retDO
end
dmul(dol::Num, dor::Num, v2d::Bijection{Num, Num}; use_asir=false) = dmul(dol, value(dor), v2d; use_asir=use_asir)
dmul(dol, dor, v2d::Bijection{Num, Num}; use_asir=false) = dol*dor