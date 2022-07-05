
# private, do not check zero-dimensionality of ideal and validity of std_mons
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

struct PfaffianSystem 
	A::Dict{Num, Matrix{Num}}
	std_mons::Vector{Num}
	v2d::Bijection{Num, Num}
end

function PfaffianSystem(I::DIdeal)
	return computePfaffSys(I, I.v2d.domain |> collect |> sort |> OrderedSet)
end