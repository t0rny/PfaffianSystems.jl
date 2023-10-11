using PfaffianSystems
using DataStructures
using Bijections
using Test

@testset "PfaffianSystems.jl" begin
    # Write your tests here.
end

# using PfaffianSystems: makeTestVarsAndIdeal

using AbstractAlgebra: QQ

@testset "DiffOps.jl" begin
    D, (x, y), (dx, dy) = weyl_algebra(["x", "y"])

    @test isequal(gens(D), [x, y])          # gens
    @test isequal(dgens(D), [dx, dy])       # gens

    @test isequal(vars(x*dx+y^2), [x, y])   # vars
    @test isequal(vars(dx+dy^2), [])        # vars with no var
    @test isequal(dvars(x*dx+y^2), [dx])    # dvars
    @test isequal(dvars(x^2+y^2), [])       # dvars with no dvar

    @test isvar(x)          # isvar
    @test !isvar(dx)        # isvar with dvar
    @test !isvar(x + dx)    # isvar with differential operator
    @test isdvar(dx)        # isdvar
    @test !isdvar(x)        # isdvar with var
    @test !isdvar(x + dx)    # isdvar with differential operator

    @test isequal(evaluate(x*y + dx*dy, [x, dx], [y, dy]), y^2 + dy^2)  # evaluate
end

@testset "WeylAlgebra.jl" begin
    @test_nowarn weyl_algebra("x")      # generate 1-d Weyl algebra
    D, x, dx = weyl_algebra("x")
    @test isequal(x*x, x^2)             # power of variable
    @test isequal(dx*x, x*dx+1)         # Leibniz rule
    @test isequal(dx*dx, dx^2)          # power of derivative
    @test isequal(x^2*x, x*x^2)         # commutativity between product and power
    @test isequal((x*dx)*x, x*(dx*x))   # associativity of product

    @test_nowarn weyl_algebra(["x","y"]) # generate 2-d Weyl algebra
    D, (x,y), (dx,dy) = weyl_algebra(["x","y"])
    @test isequal(vars(x*dx*dy), vars(x))
    @test isequal(dvars(x*y*dx), dvars(dx))
    @test isequal(vars(x*y*dx), vars(x + y))
    @test isequal(dvars(x*dx*dy), dvars(dx + dy))

    @test isequal(x*y, y*x)     # commutativity of variables
    @test isequal(dx*dy, dy*dx) # commutativity of derivatives
    @test isequal(y*dx, dx*y)   # commutativity of variables and derivatives 
    @test isequal(x*(y + dy), x*y + x*dy) # distributivity of product with variables over sum
    @test isequal(dx*(x + y), dx*x + dx*y) # distributivity of product with derivatives over sum
end

@testset "DiffOpRing.jl" begin
    @test_nowarn diff_op_ring("x")      # generate 1-d Weyl algebra
    D, x, dx = diff_op_ring("x")
    @test isequal(x*x, x^2)             # power of variable
    @test isequal(dx*x, x*dx+1)         # Leibniz rule
    @test isequal(dx*dx, dx^2)          # power of derivative
    @test isequal(x^2*x, x*x^2)         # commutativity between product and power
    @test isequal((x*dx)*x, x*(dx*x))   # associativity of product

    @test_nowarn diff_op_ring(["x","y"]) # generate 2-d Weyl algebra
    D, (x,y), (dx,dy) = diff_op_ring(["x","y"])
    @test isequal(vars(x*dx*dy), vars(x))
    @test isequal(dvars(x*y*dx), dvars(dx))
    @test isequal(vars(x*y*dx), vars(x + y))
    @test isequal(dvars(x*dx*dy), dvars(dx + dy))

    @test isequal(x*y, y*x)     # commutativity of variables
    @test isequal(dx*dy, dy*dx) # commutativity of derivatives
    @test isequal(y*dx, dx*y)   # commutativity of variables and derivatives 
    @test isequal(x*(y + dy), x*y + x*dy) # distributivity of product with variables over sum
    @test isequal(dx*(x + y), dx*x + dx*y) # distributivity of product with derivatives over sum
    @test isequal((x+dx)//y, x//y + dx//y) # distributivity of division
    @test isequal((x+dx)*y^-1, x//y + dx//y) # distributivity of division
    @test isequal(x//x, one(D)) # reduction
end

@testset "Coersion" begin
    D1, x, dx = weyl_algebra("x")
    R1, X, dX = diff_op_ring("x")
    D2, (x1, x2), (dx1, dx2) = weyl_algebra(["x", "y"])
    R2, (X1, X2), (dX1, dX2) = diff_op_ring(["x", "y"])

    @test isequal(coerce(dx*x+1, D2) |> parent, D2)
    @test isequal(coerce(dX*X+1, R2) |> parent, R2)    
    @test isequal(coerce(dx*x+1, R1) |> parent, R1)    
    @test isequal(coerce(dx1*x1+dx2*x2, R2) |> parent, R2)
end

@testset "DiffOpRings.jl" begin
    D, (x,y,z), (dx,dy,dz) = diff_op_ring(["x","y","z"])

    # r is equalt to 0
    f = dx*dy^3
    g = [x*dy+1, dx]
    r, q = normalform(f, g)
    @test isequal(f, q[1] * g[1] + q[2] * g[2] + r)

    # r is not equal to 0 and q[1] includes a non-constant coefficient
    f = dx*dy^3
    g = [x*dy+1, dx]
    r, q = normalform(f, g)
    @test isequal(f, q[1] * g[1] + q[2] * g[2] + r)

    # Example 6.1.6 in [N. Takayama, ``Chapter 6: Grobner Basis for Rings of Differnetial Operators and Applications,'', Grobner Bases, Hibi eds., Springer, 2013]
    f = dx*dy^3
    g = [dx*dy+1, 2*y*dy^2-dx+3*dy+2*x]
    r , q =  normalform(f, g)
    @test isequal(f, q[1] * g[1] + q[2] * g[2] + r)
end


