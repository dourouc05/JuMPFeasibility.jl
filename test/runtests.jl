using Test
using LinearAlgebra
using SparseArrays

using JuMP
using JuMPFeasibility

@testset "JuMPFeasibility" begin
    m = Model()
    @variable(m, x[1:2])
    @constraint(m, c1, sum(x) <= 2.0)
    @constraint(m, c2, sum(x) >= 0.0)
    @constraint(m, c3, x[1] == 0.5)
    @constraint(m, c4, 0.5 <= x[2] <= 1.0)
    @constraint(m, c5, x[1] == 1.5)

    vars = Dict(x[1] => 0.5, x[2] => 1.0)
    var_to_idx = Dict(x[1] => 1, x[2] => 2)
    f = vidx -> vars[vidx]

    ## is_satisfied
    @test is_satisfied(c1, f)
    @test is_satisfied(c2, f)
    @test is_satisfied(c3, f)
    @test is_satisfied(c4, f)
    @test !is_satisfied(c5, f)

    # Increase the tolerance.
    @test is_satisfied(c5, f, ε=1.0)

    ## is_tight
    @test !is_tight(c1, f)
    @test !is_tight(c2, f)
    @test is_tight(c3, f)
    @test !is_tight(c4, f)
    @test !is_tight(c5, f)

    # Consider the interval constraint to be degenerate by increasing the tolerance. 
    @test is_tight(c4, f, ε=1.0)

    ## coefficients
    a1 = coefficients(c1, var_to_idx)
    a2 = coefficients(c2, var_to_idx)
    a3 = coefficients(c3, var_to_idx)
    a4 = coefficients(c4, var_to_idx)
    a5 = coefficients(c5, var_to_idx)

    @test a1 == [1.0, 1.0]
    @test a2 == [1.0, 1.0]
    @test a3 == [1.0, 0.0]
    @test a4 == [0.0, 1.0]
    @test a5 == [1.0, 0.0]
end