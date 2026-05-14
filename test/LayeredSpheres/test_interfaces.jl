using Test
using MeanFieldHom
using TensND
using LinearAlgebra

# =============================================================================
#  Imperfect interfaces: spring ( 2014) and membrane
#  (Gurtin-Murdoch 1975) — bulk mode.
# =============================================================================

@testset "SpringInterface — limits recover (bulk)" begin
    κ₀, μ₀ = 100.0, 70.0
    C₀ = TensISO{3}(3κ₀, 2μ₀)
    C₁ = TensISO{3}(3κ₀ * 2, 2μ₀ * 2)

    s_perfect = LayeredSphere((0.5, 1.0), (C₁, C₀))
    α_perfect = MeanFieldHom.LayeredSpheres._bulk_localization(s_perfect, κ₀, μ₀)

    # k → 0 ≡ perfect interface (normal-only)
    intf_tight = (SpringInterface(1.0e-14), PerfectInterface{Float64}())
    s_tight = LayeredSphere((0.5, 1.0), (C₁, C₀); interfaces = intf_tight)
    α_tight = MeanFieldHom.LayeredSpheres._bulk_localization(s_tight, κ₀, μ₀)
    for k in 1:2
        @test α_tight[k] ≈ α_perfect[k] rtol = 1.0e-10
    end

    # k → ∞ : core decoupled → core bulk localization → 0
    intf_loose = (SpringInterface(1.0e9), PerfectInterface{Float64}())
    s_loose = LayeredSphere((0.5, 1.0), (C₁, C₀); interfaces = intf_loose)
    α_loose = MeanFieldHom.LayeredSpheres._bulk_localization(s_loose, κ₀, μ₀)
    @test abs(α_loose[1]) < 1.0e-6
end

@testset "SpringInterface — monotone in compliance" begin
    κ₀, μ₀ = 100.0, 70.0
    C₀ = TensISO{3}(3κ₀, 2μ₀)
    C₁ = TensISO{3}(3κ₀ * 2, 2μ₀ * 2)

    ks = [1.0e-5, 1.0e-3, 1.0e-1, 1.0]
    αs = Float64[]
    for k in ks
        s = LayeredSphere(
            (0.5, 1.0), (C₁, C₀);
            interfaces = (SpringInterface(k), PerfectInterface{Float64}())
        )
        push!(αs, MeanFieldHom.LayeredSpheres._bulk_localization(s, κ₀, μ₀)[1])
    end
    for i in 1:(length(ks) - 1)
        @test αs[i] ≥ αs[i + 1]
    end
end

@testset "SpringInterface(kn, kt) two-compliance constructor" begin
    s_norm = SpringInterface(0.01)         # convenience (kt = 0)
    s_full = SpringInterface(0.01, 0.02)   # two compliances
    @test s_norm.kn ≈ 0.01
    @test s_norm.kt == 0.0
    @test s_full.kn ≈ 0.01
    @test s_full.kt ≈ 0.02
end

@testset "MembraneInterface — Gurtin-Murdoch bulk limit" begin
    κ₀, μ₀ = 100.0, 70.0
    C₀ = TensISO{3}(3κ₀, 2μ₀)
    C₁ = TensISO{3}(3κ₀ * 2, 2μ₀ * 2)

    # κs = μs = 0 recovers perfect interface.
    s_perfect = LayeredSphere((0.5, 1.0), (C₁, C₀))
    α_perfect = MeanFieldHom.LayeredSpheres._bulk_localization(s_perfect, κ₀, μ₀)
    intf_zero = (MembraneInterface(0.0, 0.0), PerfectInterface{Float64}())
    s_zero = LayeredSphere((0.5, 1.0), (C₁, C₀); interfaces = intf_zero)
    α_zero = MeanFieldHom.LayeredSpheres._bulk_localization(s_zero, κ₀, μ₀)
    for k in 1:2
        @test α_zero[k] ≈ α_perfect[k] rtol = 1.0e-12
    end

    # Positive κs reinforces the interface: bulk localization in layer 1
    # should differ from perfect case.
    intf_m = (MembraneInterface(1.0, 0.5), PerfectInterface{Float64}())
    s_m = LayeredSphere((0.5, 1.0), (C₁, C₀); interfaces = intf_m)
    α_m = MeanFieldHom.LayeredSpheres._bulk_localization(s_m, κ₀, μ₀)
    @test α_m[1] != α_perfect[1]
end

@testset "SpringInterface eltype inference" begin
    s = LayeredSphere(
        (0.5, 1.0), (TensISO{3}(1.0, 1.0), TensISO{3}(1.0, 1.0));
        interfaces = (SpringInterface(0.01), PerfectInterface{Float64}())
    )
    @test layer_interface(s, 1) isa SpringInterface
    @test layer_interface(s, 2) isa PerfectInterface
    @test eltype(layer_interface(s, 1)) === Float64
end
