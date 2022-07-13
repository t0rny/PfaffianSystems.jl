# --------------------- Utilities ---------------------
function removeSqBra(eq::AbstractString)
	return replace(eq, r"[\[, \]]"=>"")
end

function add_ast(eq::AbstractString)
	return replace(eq, r"([0-9]+)([\(a-zA-Z])"=>s"\1*\2")
end

function addSqBra(eq::AbstractString)
	return replace(eq, r"([a-zA-Z]+)([0-9]+)"=>s"\1[\2]")
end

function slash2Dslash(eq::AbstractString)
	return replace(eq, r"([0-9a-zA-Z\)])\/([0-9a-zA-Z\(])"=>s"\1//\2")
end


"""
	vec2str(v::Vector{Num}; delim=",")
	vec2str(v::Symbolics.Arr; delim=",")
	vec2str(v::OrderedSet{Num}; delim=",")
Return a string consiting of all elements of `v` with delimiter `delim`. 
"""
vec2str(v::AbstractVector; delim=",") = string.(v) |> (s->join(s, delim))
vec2str(v::Vector{Num}; delim=",") = string.(v) .|> removeSqBra .|> add_ast |> (s->join(s, delim))
vec2str(v::Symbolics.Arr; delim=",") = vec2str(v |> scalarize; delim=delim)
vec2str(v::OrderedSet{Num}; delim=",") = vec2str(v |> collect; delim=delim)
function vec2str(vs...; delim=",")
	strs = Vector{String}()
	for v in vs
		# if typeof(v) <: Symbolics.Arr
		# 	# push!(strs, join(removeSqBra.(string.(scalarize(v))), delim))
		# 	v = v |> scalarize .|> string 
		# else
		# 	# push!(strs, join(removeSqBra.(string.(v)), delim))
		# 	v = v .|> string 
		# end
		v |> vec2str |> (s->push!(strs, s))
	end
	return join(strs, delim)
end


"""
Return the version of available Asir as an integer if exists, and otherwise return `nothing`. 
"""
function isAsirAvailable()
	asir_res = runAsir("version();") |> parseAsir
	if length(asir_res) == 2 
		asir_ver = tryparse(Int, asir_res[1])
		if !isnothing(asir_ver)
			return asir_ver
		end
	end

	# println("Worning: Asir is not available")
	return nothing
end


# --------------------- Interfaces ---------------------
"""
	runAsir(commands::AbstractString)

Run `commands` on Asir. The raw response of Asir is returned as a string. 
"""
function runAsir(commands::AbstractString; errMsg=false)
	commands = replace(commands, r"[\n\t]"=>"") |> add_ast
	# commands = replace(commands, r"([0-9]+)([a-zA-Z])"=>s"\1*\2")
	tmpFile, tmpIO = mktemp()
	tmpErrFile, tmpErrIO = mktemp()
	# tmpErrFile, tmpErrIO = mktemp()
	# @info tmpFile
	print(tmpIO, commands)
	close(tmpIO)
	close(tmpErrIO)
	# close(tmpErrIO)
	# open("./asir.tmp", "w") do out
	# 	print(out, commands)
	# end
	# o = read(`asir -quiet -f ./asir.tmp`, String)
	# return addSqBra(o)
	# return read(`asir -quiet -f $tmpFile '2>&1'`, String)

	# redirect stderror message to temporary file as it makes tests fail
	asirout = read(pipeline(`asir -quiet -f $tmpFile`, stderr=tmpErrFile), String)
	if errMsg
		@info open(tmpErrFile, "r") do f read(f, String) end
	end
	# @info open(tmpErrFile, "r") do f read(f, String) end
	return asirout
end


"""
	parseAsir(asir_res::String)

Add square brackets and separete the response into lines. 
"""
function parseAsir(asir_res::AbstractString)
	return asir_res |> (s->split(s, "\n"))
end

function evalAsir(asir_res::AbstractString, vars_list::Vector{Num})
	eval(Meta.parse("asir_tmpFunc($(vec2str(vars_list))) = $(asir_res |> slash2Dslash)"))
	# return Base.invokelatest(asir_tmpFunc, vars_list...)
	return (@invokelatest asir_tmpFunc(vars_list...)) .|> Num
end

function asir_derivative(sym::Num, var::Num)
	vars_list = get_variables(sym) .|> Num
	asir_cmd = "diff($sym, $var);"
	# """
	# diff($sym, $var);
	# """
	asir_res = asir_cmd |> runAsir |> parseAsir
	return evalAsir(asir_res[1], vars_list)
end
asir_derivative(syms::AbstractArray{Num}, var::Num) = asir_derivative.(syms, var)

function asir_reduce(sym::Num)
	vars_list = get_variables(sym) .|> Num
	asir_cmd = "red($sym);"
	asir_res = asir_cmd |> runAsir |> parseAsir
	return evalAsir(asir_res[1], vars_list)
end