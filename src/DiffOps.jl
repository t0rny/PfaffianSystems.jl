
# PolyDiffOp(a::Vector{T}, x::MonomialVector{false}) = Polynomial{false, Rational{Integer}}(a, x)

const DGen = PolyVar{false}
const DGenMon = Monomial{false}

# Use common dictionaries for maintaining relations: 
# x <--> dx
# "x" --> (x as DGen, x as Num)
const VAR2DIFF = Bijection{DGen, DGen}()
const NAME2SYM = OrderedDict{String, Pair{DGen, Num}}()


function get_var2diff()
	return VAR2DIFF
end

function get_name2sym()
	return NAME2SYM
end

get_pvar_from_name(n2s, name) = n2s[name].first
get_sym_from_name(n2s, name) = n2s[name].second

function clear_dicts()
	v2d = get_var2diff()
	n2s = get_name2sym()

	empty!(v2d)
	empty!(n2s)
	nothing
end

function check_dicts_sanity()
	v2d = get_var2diff()
	n2s = get_name2sym()

	for v in union(domain(v2d))
		# if !haskey(n2s, DP.name(v))
		# 	@assert "$(v) is not included in keys of PVAR2SYM"
		# end
		@assert haskey(n2s, DP.name(v)) "$(v) is not included in keys of NAME2SYM"
	end
	# @assert p.v2d === q.v2d "v2d mismatch"
	# @assert p.p2s === q.p2s "v2d mismatch"

	# return p.v2d, p.p2s
	return true
end

function check_var_existence(name::AbstractString)
	v2d = get_var2diff()
	return check_dicts_sanity() && (name in DP.name.(domain(v2d)))
end

abstract type DiffOp end

Base.print(io::IO, dop::DiffOp) = print(io, dop.p)
Base.show(io::IO, dop::DiffOp) = show(io, dop.p)

# struct Dgen <: DiffOp
	# p::PolyVar{false}
# 	v2d::Bijection{PolyVar{false}, PolyVar{false}}
# 	p2s::Bijection{PolyVar{false}, Num}
# end

struct PolyDiffOp <: DiffOp
	p::Polynomial{false, Rational}
# 	v2d::Bijection{DGen, DGen}
# 	p2s::Bijection{DGen, Num}

	PolyDiffOp(p::Polynomial) = new(p)
	PolyDiffOp(r::Real) = new(one(Polynomial{false, Rational})*convert(Rational, r))
	PolyDiffOp(g::DGen) = new(convert(Polynomial{false, Rational}, g))
# 	# PolyDiffOp(g::DGen, v2d, p2s) = new(convert(Polynomial{false, Rational}, g), v2d, p2s)
end


Base.copy(p::PolyDiffOp) = new(copy(p.p))
Base.isequal(p::PolyDiffOp, q::PolyDiffOp) = isequal(p.p, q.p)

Base.:+(p::PolyDiffOp, q::PolyDiffOp) = PolyDiffOp(p.p + q.p)
Base.:+(p::PolyDiffOp, r::Real) = PolyDiffOp(p.p + convert(Rational, r))
Base.:+(r::Real, q::PolyDiffOp) = PolyDiffOp(convert(Rational, r) + q.p)

Base.:-(p::PolyDiffOp, q::PolyDiffOp) = PolyDiffOp(p.p - q.p)
Base.:-(p::PolyDiffOp, r::Real) = PolyDiffOp(p.p - convert(Rational, r))
Base.:-(r::Real, q::PolyDiffOp) = PolyDiffOp(convert(Rational, r) - q.p)

Base.:*(p::PolyDiffOp, q::PolyDiffOp) = PolyDiffOp(p.p * q.p) |> canonicalize
Base.:*(p::PolyDiffOp, r::Real) = PolyDiffOp(p.p*convert(Rational, r))
Base.:*(r::Real, p::PolyDiffOp) = PolyDiffOp(p.p*convert(Rational, r))

Base.:^(p::PolyDiffOp, n::Integer) = PolyDiffOp(p.p^n) |> canonicalize

Base.one(::Type{PolyDiffOp}) = PolyDiffOp(one(Polynomial{false, Rational}))
Base.zero(::Type{PolyDiffOp}) = PolyDiffOp(zero(Polynomial{false, Rational}))

# TODO: better to implement differential operators with coefficients in rational functions 
# struct RatDiffOp <: DiffOp
# 	p::Polynomial{false, RationalPoly{Polynomial{true, Rational}, Polynomial{true, Rational}}}

# 	RatDiffOp(p::Polynomial) = new(p)
# end

# const PolyDiffOp = Polynomial{false, Rational{Int64}}
# PolyDiffOp(p::Polynomial) = convert(PolyDiffOp, p) |> canonicalize
# PolyDiffOp(g::DGen) = convert(PolyDiffOp, g)

# Refer to construction of PolyForm in [SymbolicUtils](https://github.com/JuliaSymbolics/SymbolicUtils.jl)

# let const VAR2DIFF = Bijection{DGen, DGen}(), const PV2SYM = Bijection{Dgen, Num}()
# end

# const VAR2DIFF = Ref(WeakRef(nothing))
# const PV2SYM = Ref(WeakRef(nothing))
# function get_varorder()
# 	return VARORDER
# end

# function get_var2diff()
# 	v = VAR2DIFF[].value
# 	if v === nothing
# 		d = Bijection{DGen, DGen}()
# 		VAR2DIFF[] = WeakRef(d)
# 		return d
# 	else
# 		return v
# 	end
# end

# function get_pv2sym()
# 	v = PV2SYM[].value
# 	if v === nothing
# 		d = Bijection{DGen, Num}()
# 		PV2SYM[] = WeakRef(d)
# 		return d
# 	else
# 		return v
# 	end
# end

# TODO: differential operators with coefficients in rational functions should be implemented
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

# function genVars2(name::AbstractString, v2d=get_var2diff(), p2s=get_pv2sym())
function genVars2(name::AbstractString)
	v2d = get_var2diff()
	n2s = get_name2sym()
	# vorder = get_varorder()

	@assert !startswith(name, 'd') "variable name \"$name\" must start from letters except for \"d\""
	# @assert !in(name, DP.name.(v2d.domain)) "already exists variable named \"$name\""

	# s2p = active_inv(p2s)

	# if haskey(s2p, symvar) && haskey(s2p, symdop) && v2d[s2p[symvar]] == s2p[symdop]
	if check_var_existence(name) 
		@warn "variable already exist"
		pvar = get_pvar_from_name(n2s, name)
		return PolyDiffOp(pvar), PolyDiffOp(v2d[pvar])
	else

		var_ex = Symbol(name)
		(symvar,) = @variables $var_ex

		pvar = DGen(name)
		pdop = DGen('d'*name)

		# s2p[symvar] = pvar
		# s2p[symdop] = pdop
		n2s[name] = pvar=>symvar
		v2d[pvar] = pdop
		# push!(vorder, pvar)
		return PolyDiffOp(pvar), PolyDiffOp(pdop)
	end
end
# genVars2(name::AbstractString) = genVars2(name, get_var2diff(), get_pv2sym())

function genVars2(name::AbstractString, n::Integer)
	pvars = Vector{PolyDiffOp}(undef, n)
	pdiffops = Vector{PolyDiffOp}(undef, n)

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

function (dop::DiffOp)(F::Num)
	coeffs = DP.coefficients(dop.p)
	mons = DP.monomials(dop.p)
	appliedNums = Vector{Num}(undef, length(mons))
	for i in eachindex(mons)
		appliedNums[i] = mons[i](F)
	end
	return coeffs'*appliedNums
end

function (dopmon::DGenMon)(F::Num)
	d2v = get_var2diff() |> active_inv
	n2s = get_name2sym()

	retNum = copy(F)

	ss = DP.variables(dopmon)
	es = DP.exponents(dopmon)

	# dopmon should be canonical 
	for i in reverse(eachindex(ss))
		if isDiff(ss[i])
			for _ in 1:es[i]
				retNum = derivative(retNum, get_sym_from_name(n2s, DP.name(d2v[ss[i]])))
			end
		else
			retNum = get_sym_from_name(n2s, DP.name(ss[i]))^es[i]*retNum
		end
	end
	return retNum
end

# function (g::DGen)(F, n2s, d2v)
# 	derivative(F, get_sym_from_name(n2s, d2v[g]))
# end

# function apply_do2(DiffOp::PolyDiffOp, F::Num, v2d::{})
# end

# isDiff(s::DGen) = startswith(DP.name(s), 'd')
isVar(s::DGen) = s in domain(get_var2diff())
isDiff(s::DGen) = s in image(get_var2diff())
function canonicalize(diffop::PolyDiffOp)
	coeffs = DP.coefficients(diffop.p)
	mons = DP.monomials(diffop.p)
	canOps = Vector{PolyDiffOp}(undef, length(mons))
	for i in eachindex(mons)
		canOps[i] = canonicalize(mons[i])
	end
	return coeffs'*canOps
end
function canonicalize(dopmon::DGenMon)
	# v2d = get_var2diff()
	# ordVars = get_name2sym() |> values |> collect .|> first
	# vlist = DP.variables(dopmon)
	# elist = DP.exponents(dopmon)
	# for v in ordVars
	# 	v in vlist || continue
	# 	retPoly= one(Polynomial{false, Rational})
	# 	vindcs = @. (vlist == v) || (vlist == v2d[v])
	# 	vs = vlist[vindcs]
	# 	es = elist[vindcs]
	# end
	d2v =  get_var2diff() |> active_inv
	v = DP.variables(dopmon)
	e = DP.exponents(dopmon)
	# d2v =  active_inv(v2d)
	# retPolydop = PolyDiffOp(1)
	retPoly= one(Polynomial{false, Rational})
	for i in reverse(eachindex(v))
		e[i] == 0 && continue
		if isDiff(v[i])
			for _ in 1:e[i]
				retPoly= retPoly*v[i] + DP.differentiate(retPoly, d2v[v[i]])
			end
		else
			retPoly= v[i]^e[i]*retPoly
		end
	end
	# return vars, exps
	# return PolyDiffOp(retPoly)
	return [tidyup_commutatives(m) for m in DP.monomials(retPoly)] |> sum |> PolyDiffOp
end
# need to sort commutative symbols for consistent expression
function tidyup_commutatives(dopmon::DGenMon)
	v2d = get_var2diff() 
	# d2v = get_var2diff() |> active_inv
	ordVars = get_name2sym() |> values |> collect .|> first
	syms = DP.variables(dopmon)
	exps = DP.exponents(dopmon)
	vars = Vector{DGen}()
	varexps = Vector{Int64}()
	dops = Vector{DGen}()
	dopexps = Vector{Int64}()
	sizehint!(vars, length(syms))
	sizehint!(varexps, length(syms))
	sizehint!(dops, length(syms))
	sizehint!(dopexps, length(syms))

	for v in ordVars
		vinds = findall(syms .== v)
		e = sum(exps[vinds])
		if e != 0
			push!(vars, v)
			push!(varexps, e)
		end
		dinds = findall(syms .== v2d[v])
		e = sum(exps[dinds])
		if e != 0
			push!(dops, v2d[v])
			push!(dopexps, e)
		end
	end

	return DGenMon([vars; dops], [varexps; dopexps])
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

function makeTestEnvs()
	# global PS = PfaffianSystems
	# global DP = DynamicPolynomials
	x, dx = genVars2("x", 2)
	p = (x[1] + dx[1])^3 
	q = (x[2] + dx[2])^2
	# m = DP.monomials(p.p) |> DP.canonical |> (s->s[8])
	return x, dx, p, q 
end