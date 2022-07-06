using PfaffianSystems
using DataStructures
using Bijections
using Test

@testset "PfaffianSystems.jl" begin
    # Write your tests here.
end

using PfaffianSystems: makeTestVarsAndIdeal

@testset "DiffOps.jl" begin
    x, dx, var2do = genVars("x", 3)
    @test length(x) == length(dx) == length(var2do) == 3
    @test map(s->haskey(var2do, s), x) |> all
    @test map(s->haskey(inv(var2do), s), dx) |> all

    y, dy, var2do = genVars!("y", var2do)
    @test length(var2do) == 4
    @test haskey(var2do, y)
    @test haskey(inv(var2do), dy)

    @test apply_do(dx[1], x[1], var2do) == 1
    @test apply_do(dx[2]^2 + 1, sin(x[2]), var2do) == 0
end

@testset "DIdeals.jl" begin
    x, dx, var2do = genVars("x", 3)
    I = @test_nowarn DIdeal([dx[1]^2 + 1, x[2]*dx[2] - 2, x[3]*dx[3] - 1], var2do)
    y, dy, var2do = genVars!("y", 2, var2do)
    J = DIdeal([dx[1]^2 + 1, x[2]*dy[2] - 2], var2do)
    @test isequal(stdmon!(I, OrderedSet(x)), [dx[1], 1])
    @test isZeroDimensional(I)
    @test !isZeroDimensional(J)
    I3 = eliminationIdeal(I, x[1:2])
    @test isequal(I3.gens, [x[3]*dx[3] - 1])
    @test isequal(I3.v2d, Bijection(x[3], dx[3]))
    @test_nowarn intersectionIdeal(I, J)
end

using Symbolics: @variables
@testset "AsirWrapper.jl" begin
    x, dx, var2do = genVars("x", 3)
    @variables y[1:2]
    os = OrderedSet((x[2], x[1]))

    @test !isnothing(isAsirAvailable())
    @test vec2str(x) == "x1,x2,x3"
    @test vec2str(y) == "y1,y2"
    @test vec2str(os) == "x2,x1"
end

@testset "PfaffSys.jl" begin
    x, dx, var2do, I = makeTestVarsAndIdeal()
    @test_nowarn PfaffianSystem(I)
end