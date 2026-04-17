using Test
using MeanFieldHom

@testset "MeanFieldHom" begin
    @testset "Core" begin
        include("Core/test_traits.jl")
        include("Core/test_newton.jl")
    end

    @testset "Elasticity" begin
        include("Elasticity/test_hill.jl")
    end

    @testset "Cracks" begin
        include("Cracks/test_cod.jl")
    end

    @testset "Conductivity" begin
        include("Conductivity/test_hill_order2.jl")
    end

    @testset "Regression" begin
        include("regression/test_hill_cases.jl")
        include("regression/test_crack_cases.jl")
    end
end
