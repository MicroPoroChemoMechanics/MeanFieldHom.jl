using Test
using MeanFieldHom

@testset "MeanFieldHom" begin
    @testset "Elliptic" begin
        include("Elliptic/test_elliptic.jl")
    end

    @testset "Core" begin
        include("Core/test_traits.jl")
        include("Core/test_newton.jl")
        include("Core/test_newton_cylinder.jl")
    end

    @testset "Elasticity" begin
        include("Elasticity/test_hill.jl")
        include("Elasticity/test_hill_cylinder.jl")
        include("Elasticity/test_shape_tensor.jl")
        include("Elasticity/test_eshelby.jl")
        include("Elasticity/test_localization.jl")
        include("Elasticity/test_contribution.jl")
        include("Elasticity/test_hill_nestedquadgk_oblate.jl")
        include("Elasticity/test_hill_ti_coaxial.jl")
    end

    @testset "Cracks" begin
        include("Cracks/test_cod.jl")
        include("Cracks/test_cod_ti_aligned.jl")
        include("Cracks/test_residue_accuracy.jl")
        include("Cracks/test_thermal.jl")
    end

    @testset "Conductivity" begin
        include("Conductivity/test_hill_order2.jl")
        include("Conductivity/test_hill_cylinder.jl")
        include("Conductivity/test_eshelby.jl")
        include("Conductivity/test_localization.jl")
    end

    @testset "LayeredSpheres" begin
        include("LayeredSpheres/test_bulk.jl")
        include("LayeredSpheres/test_interfaces.jl")
        include("LayeredSpheres/test_incompressible.jl")
        include("LayeredSpheres/test_conductivity.jl")
        include("LayeredSpheres/test_christensen.jl")
        include("LayeredSpheres/test_generic.jl")
    end

    @testset "Regression" begin
        include("regression/test_hill_cases.jl")
        include("regression/test_crack_cases.jl")
        include("regression/test_anisotropic.jl")
    end
end
