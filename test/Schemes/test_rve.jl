# =============================================================================
#  test_rve.jl
#
#  Validates the RVE / Phase / Amount data model:
#   1. Construction and progressive registration of matrix + inclusion phases.
#   2. Matrix volume fraction is implicit and ignores crack densities.
#   3. Distribution-shape coercion (nothing → UniformDistribution sphere ;
#      AbstractInclusion auto-wrapped ; passthrough for explicit
#      AbstractDistributionShape).
#   4. Argument validation: duplicate names, missing matrix, mutually
#      exclusive fraction/density kwargs, negative amounts.
#   5. Element-type propagation : Float64 (default), ForwardDiff.Dual,
#      Complex{Float64}.
# =============================================================================

using Test
using MeanFieldHom
using TensND
using ForwardDiff

@testset "RVE — basic construction & accessors" begin
    rve = RVE(:M)
    @test rve isa RVE{Float64}
    @test isempty(rve.phase_names)
    @test rve.distribution_shape isa UniformDistribution

    # Matrix
    C₀ = TensISO{3}(30.0, 10.0)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:C => C₀))
    @test matrix_phase(rve) isa Phase
    @test matrix_property(rve, :C) === C₀

    # Inclusion (volume fraction)
    C₁ = TensISO{3}(60.0, 20.0)
    add_phase!(rve, :I, Ellipsoid(1.0, 1.0, 0.3), Dict(:C => C₁); fraction = 0.2)
    @test rve.amounts[:I] isa VolumeFraction
    @test volume_fraction(rve, :I) ≈ 0.2
    @test crack_density(rve, :I) == 0.0

    # Crack (density)
    add_phase!(rve, :CRACK, PennyCrack(1.0), Dict(:C => C₀); density = 0.05)
    @test rve.amounts[:CRACK] isa CrackDensity
    @test crack_density(rve, :CRACK) ≈ 0.05
    @test volume_fraction(rve, :CRACK) == 0.0

    # Implicit matrix fraction = 1 - 0.2 (cracks excluded from sum)
    @test matrix_volume_fraction(rve) ≈ 0.8
    @test volume_fraction(rve, :M) ≈ 0.8

    # Insertion order respected
    @test inclusion_phase_names(rve) == [:I, :CRACK]
    @test rve.phase_names == [:M, :I, :CRACK]

    # No-op validation succeeds
    @test validate_rve(rve) === rve
end

@testset "RVE — argument errors" begin
    rve = RVE(:M)
    C₀ = TensISO{3}(30.0, 10.0)

    # Cannot add a phase before the matrix...
    @test_throws ArgumentError validate_rve(rve)

    add_matrix!(rve, Ellipsoid(1.0), Dict(:C => C₀))
    # Cannot add the matrix twice
    @test_throws ArgumentError add_matrix!(rve, Ellipsoid(1.0), Dict(:C => C₀))
    # Cannot use the matrix name for a phase
    @test_throws ArgumentError add_phase!(rve, :M, Ellipsoid(1.0), Dict(:C => C₀); fraction = 0.1)

    add_phase!(rve, :I, Ellipsoid(1.0), Dict(:C => C₀); fraction = 0.2)
    # Duplicate phase name
    @test_throws ArgumentError add_phase!(rve, :I, Ellipsoid(1.0), Dict(:C => C₀); fraction = 0.1)
    # Must specify exactly one of fraction / density
    @test_throws ArgumentError add_phase!(rve, :J, Ellipsoid(1.0), Dict(:C => C₀))
    @test_throws ArgumentError add_phase!(
        rve, :J, Ellipsoid(1.0), Dict(:C => C₀);
        fraction = 0.1, density = 0.1
    )
    # Negative amount caught by validate_rve
    add_phase!(rve, :J, Ellipsoid(1.0), Dict(:C => C₀); fraction = -0.1)
    @test_throws ArgumentError validate_rve(rve)
end

@testset "RVE — distribution_shape coercion" begin
    # Default → UniformDistribution(Ellipsoid(1.0))
    rve_default = RVE(:M)
    @test rve_default.distribution_shape isa UniformDistribution
    @test rve_default.distribution_shape.shape isa Ellipsoid

    # AbstractInclusion auto-wrapped
    ell = Ellipsoid(2.0, 1.0, 1.0)
    rve_wrap = RVE(:M; distribution_shape = ell)
    @test rve_wrap.distribution_shape isa UniformDistribution
    @test rve_wrap.distribution_shape.shape === ell

    # Explicit AbstractDistributionShape passthrough
    ds = UniformDistribution(Ellipsoid(1.0, 1.0, 0.5))
    rve_explicit = RVE(:M; distribution_shape = ds)
    @test rve_explicit.distribution_shape === ds
end

@testset "RVE — element-type propagation" begin
    # ForwardDiff.Dual amounts (sensitivity analysis on volume fraction)
    DT = ForwardDiff.Dual{Nothing, Float64, 1}
    rve_dual = RVE(:M; T = DT)
    @test rve_dual isa RVE{DT}
    C₀ = TensISO{3}(30.0, 10.0)
    add_matrix!(rve_dual, Ellipsoid(1.0), Dict(:C => C₀))
    add_phase!(
        rve_dual, :I, Ellipsoid(1.0), Dict(:C => C₀);
        fraction = ForwardDiff.Dual{Nothing}(0.2, 1.0)
    )
    @test rve_dual.amounts[:I] isa VolumeFraction{DT}
    @test ForwardDiff.value(matrix_volume_fraction(rve_dual)) ≈ 0.8
    @test ForwardDiff.partials(matrix_volume_fraction(rve_dual))[1] ≈ -1.0

    # Complex{Float64} amounts (frequency-domain volume-fraction sweep —
    # rarely useful, but compatibility must hold).
    rve_c = RVE(:M; T = Complex{Float64})
    add_matrix!(rve_c, Ellipsoid(1.0), Dict(:C => C₀))
    add_phase!(
        rve_c, :I, Ellipsoid(1.0), Dict(:C => C₀);
        fraction = 0.2 + 0.0im
    )
    @test rve_c.amounts[:I] isa VolumeFraction{Complex{Float64}}
    @test matrix_volume_fraction(rve_c) ≈ 0.8 + 0.0im
end

@testset "RVE — Amounts unit tests" begin
    # Construction
    @test VolumeFraction(0.3) isa VolumeFraction{Float64}
    @test CrackDensity(0.05) isa CrackDensity{Float64}

    # _sums_to_unit dispatch
    @test MeanFieldHom.Schemes._sums_to_unit(VolumeFraction(0.3))
    @test !MeanFieldHom.Schemes._sums_to_unit(CrackDensity(0.05))

    # eltype
    @test eltype(VolumeFraction(0.3)) === Float64
    @test eltype(CrackDensity{Complex{Float64}}(0.0 + 0.0im)) === Complex{Float64}

    # amount_value
    @test MeanFieldHom.Schemes.amount_value(VolumeFraction(0.3)) ≈ 0.3
    @test MeanFieldHom.Schemes.amount_value(CrackDensity(0.05)) ≈ 0.05
end
