const ASIR_DEBUG = false
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

function Dslash2slash(eq::AbstractString)
	return replace(eq, r"([0-9a-zA-Z\)])\/\/([0-9a-zA-Z\(])"=>s"\1/\2")
end


"""
	vec2str(v::AbstractVector; delim=",")
	vec2str(v::AbstractDiffOp; delim=",")
	vec2str(vs...; delim=",")
Return a string consiting of all elements of `v` with delimiter `delim`. 
"""
vec2str(v::AbstractVector; delim=",") = string.(v) |> (s->join(s, delim))
vec2str(v::AbstractDiffOp; delim=",") = string(v)
function vec2str(vs...; delim=",")
	strs = Vector{String}()
	for v in vs
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

	println("Worning: Asir is not available")
	return nothing
end


# --------------------- Interfaces ---------------------
"""
	runAsir(commands::AbstractString)

Run `commands` on Asir. The raw response of Asir is returned as a string. 
"""
function runAsir(commands::AbstractString; errMsg=ASIR_DEBUG)
	commands = replace(commands, r"[\n\t]"=>"") |> add_ast |> Dslash2slash
	if errMsg
		@info commands
	end
	tmpFile, tmpIO = mktemp()
	tmpErrFile, tmpErrIO = mktemp()
	print(tmpIO, commands)
	close(tmpIO)
	close(tmpErrIO)

	# redirect stderror message to temporary file as it makes tests fail
	asirout = read(pipeline(`asir -quiet -f $tmpFile`, stderr=tmpErrFile), String)
	if errMsg
		@info open(tmpErrFile, "r") do f read(f, String) end
	end
	return asirout
end


"""
	parseAsir(asir_res::String)

Add square brackets and separete the response into lines. 
"""
function parseAsir(asir_res::AbstractString)
	return asir_res |> (s->split(s, "\n")) .|> slash2Dslash
end

function evalAsir(asir_res::AbstractString, vars_list::Vector{<:AbstractDiffOp})
	tmpExpr = Meta.parse("asir_tmpFunc($(vec2str(vars_list))) = $(asir_res)")
	try
		eval(tmpExpr)
	catch
		@assert false "Error: somthing wrong in evaluating $(asir_res)"
	end
	return (@invokelatest asir_tmpFunc(vars_list...))
end
