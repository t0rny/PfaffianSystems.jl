
# private, do not check zero-dimensionality
function computePfaffSys(I::DIdeal, ordered_vars::OrderedSet{Num})

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
	std_mons = evalAsir(asir_res[1], vars_list)

	asir_res = asir_res[2:end]
	m = length(std_mons)
	n = length(ordered_vars)

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
	return computePfaffSys(I, I.v2d.domain |> collect |> sort |> OrderedSet)
end

function buildFuncA(pf::PfaffianSystem)
	vars = pf.v2d.domain |> collect |> sort
	[build_function(pf.A[v], vars)[1] |> eval for v in vars]
end

function integrate(pf::PfaffianSystem, init_vecs::Matrix{<:Real}, z_init::Vector{<:Real}, z_term::Vector{<:Real})
	d = length(pf.std_mons)
	@assert size(init_vecs)[1] == d "Error: invalid length of initial vectors"
	@assert length(pf.v2d.domain) == length(z_init) == length(z_term) "Error: invalid lengths of initial and terminal z vectors"

	z_traj(s) = ((z_term-z_init)*s + z_init, (z_term-z_init))

	funcA(s) = map(a->Base.invokelatest(a, s), buildFuncA(pf))

	PfODE = (q, param, s)->begin
		zs, dzds = (param)(s)
		reduce(+, funcA(zs).*dzds)*q
	end

	pf_ode = ODEProblem(PfODE, init_vecs, [0, 1], z_traj)
	sol = solve(pf_ode, abstol=1e-6, reltol=1e-6)
	return sol
end