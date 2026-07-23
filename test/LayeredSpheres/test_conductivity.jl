using Test
using MeanFieldHom
using TensND
using LinearAlgebra

# =============================================================================
#  Conductivity LayeredSphere — Y₁-harmonic gradient-gradient localization.
# =============================================================================

@testset "Conductivity — single-layer ≡ Ellipsoid Eshelby" begin
    for (k₀, k₁) in ((2.0, 5.0), (1.0, 0.3), (10.0, 10.0))
        K₀ = TensISO{3}(k₀)
        K₁ = TensISO{3}(k₁)
        s = LayeredSphere((1.0,), (K₁,))
        A_layered = gradient_gradient_loc(s, K₀; layer = 1)
        A_ell = gradient_gradient_loc(Ellipsoid(1.0), K₁, K₀)
        @test A_layered[1, 1] ≈ A_ell[1, 1] rtol = 1.0e-12
        @test A_layered[1, 1] ≈ 3 * k₀ / (2 * k₀ + k₁) rtol = 1.0e-12
    end
end

@testset "Conductivity — core-shell bulk localization sanity" begin
    k₀ = 2.0
    k_core = 5.0
    k_shell = 3.0
    K₀ = TensISO{3}(k₀)
    K_core = TensISO{3}(k_core)
    K_shell = TensISO{3}(k_shell)
    s = LayeredSphere((0.5, 1.0), (K_core, K_shell))

    α = MeanFieldHom.LayeredSpheres._cond_localization(s, k₀)
    @test all(isfinite, α)
    @test all(α .> 0)

    # Sum rule: Σ f_k · k_k · α_k should equal k_eff of the composite.
    k_eff = MeanFieldHom.LayeredSpheres._effective_conductivity(s, k₀)
    @test isfinite(k_eff) && k_eff > 0
end

@testset "Conductivity — impermeable (k=0) core is finite" begin
    # A solid, impermeable aggregate (k_core = 0) coated by a conductive shell:
    # the transfer-matrix `q̂n/k` term is 0/0 for the core, but regularity
    # (B = 0) gives a finite localization.  Regression for the NaN bug.
    k₀, k_shell = 1.0, 50.0
    K₀ = TensISO{3}(k₀)
    s = LayeredSphere((1.0, 1.01), (TensISO{3}(0.0), TensISO{3}(k_shell)))

    α = MeanFieldHom.LayeredSpheres._cond_localization(s, k₀)
    @test all(isfinite, α)
    N = conductivity_contribution(s, K₀)
    @test isfinite(N[1, 1])

    # ITZ mortar (impermeable aggregate + conductive shell, MT) — matches the
    # Echoes reference value D_eff/D_cp = 0.9979 at f = 0.2, D_itz/D_cp = 50
    # (here the shell conductivity k_shell = 50).
    f_inc = 0.2 * 1.01^3
    rve = RVE(:CEMENT)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:K => TensISO{3}(1.0)))
    add_phase!(rve, :AGG, s, Dict(:K => TensISO{3}(1.0)); fraction = f_inc)
    D_eff = tr(Array(homogenize(rve, MoriTanaka(), :K))) / 3
    @test D_eff ≈ 0.9979 rtol = 1.0e-3
end

@testset "KapitzaInterface — limits" begin
    K₀ = TensISO{3}(2.0)
    K₁ = TensISO{3}(5.0)

    s_perfect = LayeredSphere((0.5, 1.0), (K₁, K₀))
    α_perfect = MeanFieldHom.LayeredSpheres._cond_localization(s_perfect, 2.0)

    # ρ → 0 recovers perfect
    intf = (KapitzaInterface(1.0e-14), PerfectInterface{Float64}())
    s = LayeredSphere((0.5, 1.0), (K₁, K₀); interfaces = intf)
    α = MeanFieldHom.LayeredSpheres._cond_localization(s, 2.0)
    for k in 1:2
        @test α[k] ≈ α_perfect[k] rtol = 1.0e-9
    end

    # ρ → ∞ decouples core (α_1 → 0)
    intf_loose = (KapitzaInterface(1.0e9), PerfectInterface{Float64}())
    s_loose = LayeredSphere((0.5, 1.0), (K₁, K₀); interfaces = intf_loose)
    α_loose = MeanFieldHom.LayeredSpheres._cond_localization(s_loose, 2.0)
    @test abs(α_loose[1]) < 1.0e-6
end

@testset "SurfaceConductiveInterface — limits" begin
    K₀ = TensISO{3}(2.0)
    K₁ = TensISO{3}(5.0)

    s_perfect = LayeredSphere((0.5, 1.0), (K₁, K₀))
    α_perfect = MeanFieldHom.LayeredSpheres._cond_localization(s_perfect, 2.0)

    # conductance = 0 recovers perfect
    intf = (SurfaceConductiveInterface(0.0), PerfectInterface{Float64}())
    s = LayeredSphere((0.5, 1.0), (K₁, K₀); interfaces = intf)
    α = MeanFieldHom.LayeredSpheres._cond_localization(s, 2.0)
    for k in 1:2
        @test α[k] ≈ α_perfect[k] rtol = 1.0e-12
    end

    # Large conductance modifies localization
    intf_large = (SurfaceConductiveInterface(1.0), PerfectInterface{Float64}())
    s_large = LayeredSphere((0.5, 1.0), (K₁, K₀); interfaces = intf_large)
    α_large = MeanFieldHom.LayeredSpheres._cond_localization(s_large, 2.0)
    @test α_large[1] != α_perfect[1]
end

@testset "conductivity_contribution — single layer matches Ellipsoid" begin
    K₀ = TensISO{3}(2.0)
    K₁ = TensISO{3}(5.0)
    s = LayeredSphere((1.0,), (K₁,))
    N_sphere = conductivity_contribution(s, K₀)
    N_ell = conductivity_contribution(Ellipsoid(1.0), K₁, K₀)
    @test N_sphere[1, 1] ≈ N_ell[1, 1] rtol = 1.0e-10
end
