using Test
using MeanFieldHom
using Random

# Load DECUHR + Integrals so the `MeanFieldHomDECUHRExt` extension activates:
# several tests cross-validate the `method = :decuhr` path against the residue
# and nested-QuadGK backends. (DECUHR is a weak dependency of MeanFieldHom.)
import DECUHR, Integrals

# Load NonlinearSolve so the `MeanFieldHomNonlinearSolveExt` extension
# activates: `test_self_consistent_nls.jl` exercises SC / ASC through
# SciML algorithms (NewtonRaphson, TrustRegion) and `AutoNonlinear`.
# (NonlinearSolve is a weak dependency of MeanFieldHom.)
import NonlinearSolve

# Load Ferrite / FerriteGmsh / Gmsh so the `MeanFieldHomFerriteExt` extension
# activates: `test_ferrite_crack.jl` exercises the finite-element crack.  They
# are weak dependencies, so a missing stack skips those tests rather than
# failing the suite.
const HAS_FERRITE = try
    @eval import Ferrite, FerriteGmsh, Gmsh
    true
catch
    false
end

# Several test files draw random operators (`test_ti_alv.jl`, `test_ortho_alv.jl`,
# `test_volterra_inverse.jl`, …).  Seed once here so a CI failure is always
# reproducible locally instead of depending on the draw.
Random.seed!(20260723)

@testset "MeanFieldHom" begin
    @testset "Elliptic" begin
        include("Elliptic/test_elliptic.jl")
    end

    @testset "Core" begin
        include("Core/test_traits.jl")
        include("Core/test_rotational_average.jl")
        include("Core/test_newton.jl")
        include("Core/test_newton_cylinder.jl")
        include("Core/test_green_dipole.jl")
    end

    @testset "Elasticity" begin
        include("Elasticity/test_hill.jl")
        include("Elasticity/test_hill_2d.jl")
        include("Elasticity/test_hill_cylinder.jl")
        include("Elasticity/test_shape_tensor.jl")
        include("Elasticity/test_eshelby.jl")
        include("Elasticity/test_localization.jl")
        include("Elasticity/test_contribution.jl")
        include("Elasticity/test_hill_nestedquadgk_oblate.jl")
        include("Elasticity/test_hill_ti_coaxial.jl")
        include("Elasticity/test_param_conversions.jl")
    end

    @testset "Cracks" begin
        include("Cracks/test_cod.jl")
        include("Cracks/test_cod_ti_aligned.jl")
        include("Cracks/test_residue_accuracy.jl")
        include("Cracks/test_H_oracle.jl")
        include("Cracks/test_thermal.jl")
        include("Cracks/test_interface_stiffness.jl")
    end

    @testset "Conductivity" begin
        include("Conductivity/test_hill_order2.jl")
        include("Conductivity/test_hill_cylinder.jl")
        include("Conductivity/test_eshelby.jl")
        include("Conductivity/test_localization.jl")
    end

    @testset "Schemes" begin
        include("Schemes/test_rve.jl")
        include("Schemes/test_dispatch.jl")
        include("Schemes/test_voigt_reuss.jl")
        include("Schemes/test_one_shot.jl")
        include("Schemes/test_maxwell_pcw.jl")
        include("Schemes/test_self_consistent.jl")
        include("Schemes/test_self_consistent_nls.jl")
        include("Schemes/test_differential.jl")
        include("Schemes/test_complex_moduli.jl")
        include("Schemes/test_dual_compat.jl")
        include("Schemes/test_parameters.jl")
        include("Schemes/test_sensitivities.jl")
        include("Schemes/test_symmetrize.jl")
        include("Schemes/test_orientation.jl")
        include("Schemes/test_loc_bundles.jl")
    end

    # The user-defined-inclusion contract spans Core (the abstractions) and
    # Schemes (every consumer), so it runs after both — the scheme kernels are
    # already compiled by then, which keeps this testset cheap.
    @testset "CustomInclusions" begin
        include("Core/test_custom_inclusion.jl")
    end

    # Finite-element inclusions: skipped when the Ferrite stack is unavailable
    # (it is a weak dependency), and slow when it is — each case meshes a ball
    # and factorizes a ~10⁵-dof system.
    if HAS_FERRITE
        @testset "FiniteElementInclusions" begin
            include("ext/test_ferrite_crack.jl")
        end
    else
        @info "Ferrite / FerriteGmsh / Gmsh unavailable — skipping the " *
            "finite-element inclusion tests."
    end

    @testset "LayeredSpheres" begin
        include("LayeredSpheres/test_bulk.jl")
        include("LayeredSpheres/test_interfaces.jl")
        include("LayeredSpheres/test_incompressible.jl")
        include("LayeredSpheres/test_conductivity.jl")
        include("LayeredSpheres/test_christensen.jl")
        include("LayeredSpheres/test_generic.jl")
        include("LayeredSpheres/test_scheme_integration.jl")
    end

    @testset "LayeredSpheroids" begin
        include("LayeredSpheroids/test_legendre.jl")
        include("LayeredSpheroids/test_coupling.jl")
        include("LayeredSpheroids/test_conductivity.jl")
        include("LayeredSpheroids/test_scheme_integration.jl")
        include("LayeredSpheroids/test_local_fields.jl")
    end

    @testset "Viscoelasticity" begin
        include("Viscoelasticity/test_symmetrize_alv.jl")
        include("Viscoelasticity/test_visco_law.jl")
        include("Viscoelasticity/test_trapezoidal.jl")
        include("Viscoelasticity/test_volterra_inverse.jl")
        include("Viscoelasticity/test_hill_alv_iso.jl")
        include("Viscoelasticity/test_schemes_alv.jl")
        include("Viscoelasticity/test_sc_alv.jl")
        include("Viscoelasticity/test_sc_alv_newton.jl")
        include("Viscoelasticity/test_layered_alv.jl")
        include("Viscoelasticity/test_ti_alv.jl")
        include("Viscoelasticity/test_ortho_alv.jl")
        include("Viscoelasticity/test_ortho_dispatch_alv.jl")
        include("Viscoelasticity/test_alv_kernel_types.jl")
        include("Viscoelasticity/test_sensitivities_alv.jl")
        include("Viscoelasticity/test_order2_alv.jl")
        include("Viscoelasticity/test_extra_schemes_alv.jl")
        include("Viscoelasticity/test_crack_schemes_alv.jl")
    end

    @testset "Regression" begin
        include("regression/test_hill_cases.jl")
        include("regression/test_crack_cases.jl")
        include("regression/test_anisotropic.jl")
    end
end
