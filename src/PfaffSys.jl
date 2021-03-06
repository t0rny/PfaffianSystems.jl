
# private, do not check zero-dimensionality
function _computePfaffSys(I::DIdeal, ordered_vars::OrderedSet{Num})

	@assert isequal(ordered_vars, I.v2d.domain) "Error: unknown variable in second argument"

	asir_cmd = 
	"""
	load("yang.rr")\$
	V = [$(vec2str(ordered_vars))]\$
	yang.define_ring(["partial", V])\$
	Gb = yang.gr([$(vec2str(I.gens))])\$
	Sm = yang.stdmon(Gb);
	Pf = yang.pf(Sm, Gb)\$
	for (I=0; I<length(V); I++) {
		print(Pf[I]);
	}
	"""

	asir_res = asir_cmd |> runAsir |> parseAsir
	asir_res = filter!((s)->(startswith(s, "[")), asir_res)

	(length(asir_res) == 0) && return nothing

	# vars_list = cat(t, vars, dt, diffops; dims=1)
	vars_list = union(I.v2d.domain, I.v2d.range) |> collect
	std_mons = evalAsir(asir_res[1], vars_list) |> reverse

	asir_res = asir_res[2:end]
	m = length(std_mons)

	A = Dict{Num, Matrix{Num}}()
	# for i = 1:n
	for (i, v) in enumerate(ordered_vars)
		A[v] = asir_res[m*(i-1)+1:m*i] |> reverse .|> (s->evalAsir(s, vars_list)) .|> reverse |> (s->vcat(s...))
	end

	return PfaffianSystem(A, std_mons, copy(I.v2d))
	# eval(Meta.parse("intersect_tmpFunc($(vec2str(vars_list))) = $(asir_res[1])"))
	# elimOrdGens = Base.invokelatest(intersect_tmpFunc, vars_list...)
end

"""
Pfaffian system
"""
struct PfaffianSystem 
	A::Dict{Num, Matrix{Num}}
	std_mons::Vector{Num}
	v2d::Bijection{Num, Num}
end

function PfaffianSystem(I::DIdeal)
	@assert isZeroDimensional(I) "Error: given ideal is not zero-dimensional"
	return _computePfaffSys(I, I.v2d.domain |> collect |> sort |> OrderedSet)
end

function denomLCM(pf::PfaffianSystem)
	tmpD = Num(1)
	for v in pf.v2d.domain
		map(pf.A[v] .|> value) do r
			if isdiv(r)
				tmpD *= r.den
			end
		end
	end
	DenomFctrs = asir_fctr(tmpD)
	reduce(*, [fctr.first for fctr in DenomFctrs])
end

function applyStdMons(pf::PfaffianSystem, F::Num; use_asir=false)
	return map(pf.std_mons) do dmon
		apply_do(dmon, F, pf.v2d; use_asir=use_asir)
	end
end

function buildFuncA(pf::PfaffianSystem)
	vars = pf.v2d.domain |> collect |> sort
	return [build_function(pf.A[v], vars)[1] |> eval for v in vars], vars
end

# private
function _integrate_core(funcA, init_vecs::AbstractMatrix{Float64}, z_init::AbstractVector{Float64}, z_term::AbstractVector{Float64})
	z_traj(s) = ((z_term-z_init)*s + z_init, (z_term-z_init))
	PfODE = (dq, q, param, s)->begin
		zs, dzds = (param)(s)
		dq .= reduce(+, funcA(zs).*dzds)*q
	end

	pf_ode = ODEProblem(PfODE, init_vecs, [0, 1], z_traj)
	sol = solve(pf_ode, abstol=1e-6, reltol=1e-6)
	@assert sol.retcode == :Success "Error: ode solver failed"
	return sol.u
end

function integrate(pf::PfaffianSystem, init_vecs::AbstractMatrix{<:Real}, z_init::Dict{Num, <:Real}, z_term::Dict{Num, <:Real})
	d = length(pf.std_mons)
	@assert size(init_vecs)[1] == d "Error: invalid length of initial vectors"
	# @assert length(pf.v2d.domain) == length(z_init) == length(z_term) "Error: invalid lengths of initial and terminal z vectors"
	@assert issetequal(pf.v2d.domain, keys(z_init)) "Error: invalid variables in initial z"
	@assert issetequal(pf.v2d.domain, keys(z_term)) "Error: invalid variables in terminal z"

	exprFuncA, vars = buildFuncA(pf)
	funcA(s) = map(exprFuncA) do fA
		@invokelatest fA(s)
	end
	zvec_init = substitute(vars, (z_init),) .|> value
	zvec_term = substitute(vars, (z_term),) .|> value

	sol = _integrate_core(
		funcA, 
		convert(Matrix{Float64}, init_vecs), 
		convert(Vector{Float64}, zvec_init), 
		convert(Vector{Float64}, zvec_term)
		)
	return sol
end

function integrate(pf::PfaffianSystem, init_vecs::AbstractMatrix{<:Real}, z_traj::Dict{Num, <:AbstractVector{<:Real}})
	d = length(pf.std_mons)
	# N = size(z_traj)[2]
	@assert size(init_vecs)[1] == d "Error: invalid length of initial vectors"
	@assert issetequal(pf.v2d.domain, keys(z_traj)) "Error: invalid variables in trajectory of z"
	# @assert reduce(==, values(z_traj) .|> length) "Error: lengths of series $(values(z_traj) .|> length) are different"
	@assert (values(z_traj) .|> length) |> (s->all(t->t==s[1], s)) "Error: lengths of series $(values(z_traj) .|> length) are different"
	N = values(z_traj) |> first |> length

	exprFuncA, vars = buildFuncA(pf)
	funcA(s) = map(exprFuncA) do fA
		@invokelatest fA(s)
	end

	vecs = fill(convert(Matrix{Float64}, init_vecs), N) 
	zvec_traj = convert(Matrix{Float64}, hcat([z_traj[v] for v in vars]...)' |> collect)

	for i = 1:N-1
		zvec_init = @view zvec_traj[:, i]
		zvec_term = @view zvec_traj[:, i+1]
		sol = _integrate_core(
			funcA, 
			vecs[i], 
			zvec_init, 
			zvec_term
			)
		vecs[i+1] = sol[length(sol)]
	end
	return vecs
end

# deprecated
function integrate(pf::PfaffianSystem, init_vecs::Matrix{<:Real}, z_init::Vector{<:Real}, z_term::Vector{<:Real})
	d = length(pf.std_mons)
	@assert size(init_vecs)[1] == d "Error: invalid length of initial vectors"
	@assert length(pf.v2d.domain) == length(z_init) == length(z_term) "Error: invalid lengths of initial and terminal z vectors"

	# z_traj(s) = ((z_term-z_init)*s + z_init, (z_term-z_init))

	# funcA(s) = map(a->Base.invokelatest(a, s), buildFuncA(pf))
	exprFuncA, vars = buildFuncA(pf)
	funcA(s) = map(exprFuncA) do fA
		@invokelatest fA(s)
	end

	# PfODE = (dq, q, param, s)->begin
	# 	zs, dzds = (param)(s)
	# 	dq .= reduce(+, funcA(zs).*dzds)*q
	# end

	# pf_ode = ODEProblem(PfODE, init_vecs, [0, 1], z_traj)
	# sol = solve(pf_ode, abstol=1e-6, reltol=1e-6)
	# return sol
	sol = _integrate_core(
		funcA, 
		convert(Matrix{Float64}, init_vecs), 
		convert(Vector{Float64}, z_init), 
		convert(Vector{Float64}, z_term)
		)
	return sol
end

# deprecated
function integrate(pf::PfaffianSystem, init_vecs::Matrix{<:Real}, z_traj::Matrix{<:Real})
	d = length(pf.std_mons)
	N = size(z_traj)[2]
	@assert size(init_vecs)[1] == d "Error: invalid length of initial vectors"
	@assert length(pf.v2d.domain) == size(z_traj)[1] "Error: invalid lengths of z vectors"

	vecs = fill(convert(Matrix{Float64}, init_vecs), N) 
	z_traj = convert(Matrix{Float64}, z_traj)

	# funcA(s) = map(a->Base.invokelatest(a, s), buildFuncA(pf))
	funcA(s) = map(buildFuncA(pf)) do fA
		@invokelatest fA(s)
	end

	for i = 1:N-1
		z_init = @view z_traj[:, i]
		z_term = @view z_traj[:, i+1]
		sol = _integrate_core(
			funcA, 
			vecs[i], 
			z_init, 
			z_term
			)
		vecs[i+1] = sol[length(sol)]
	end
	return vecs
end