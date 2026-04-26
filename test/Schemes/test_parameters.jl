# =============================================================================
#  test_parameters.jl — Lentilles paramétriques (get_param / set_param).
#
#  Vérifie le round-trip get_param ∘ set_param sur chaque sous-type
#  d'AbstractParameter, la promotion correcte du type d'élément du RVE,
#  et l'absence de mutation de l'instance d'origine.
# =============================================================================

using Test
using MeanFieldHom
using TensND
using ForwardDiff

@testset "AmountParameter — round-trip & type promotion" begin
    rve = RVE(:M)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:C => TensISO{3}(30.0, 10.0)))
    add_phase!(rve, :I, Ellipsoid(1.0), Dict(:C => TensISO{3}(60.0, 20.0));
               fraction = 0.2)
    add_phase!(rve, :C, Ellipsoid(1.0, 1.0, 0.0), Dict(:C => TensISO{3}(60.0, 20.0));
               density = 0.1)

    p_f = amount(:I)
    p_e = amount(:C)
    p_M = amount(:M)

    # get_param
    @test get_param(rve, p_f) ≈ 0.2
    @test get_param(rve, p_e) ≈ 0.1
    @test get_param(rve, p_M) ≈ 1.0 - 0.2     # matrix is 1 - Σ f_inc

    # set_param round-trip (Float64)
    rve2 = set_param(rve, p_f, 0.3)
    @test get_param(rve2, p_f) ≈ 0.3
    @test eltype(rve2) === Float64

    # set_param promotes amount eltype to Dual when value is Dual
    dx = ForwardDiff.Dual{Nothing, Float64, 1}(0.2, ForwardDiff.Partials((1.0,)))
    rve3 = set_param(rve, p_f, dx)
    @test eltype(rve3) <: ForwardDiff.Dual

    # CrackDensity is preserved as CrackDensity, not VolumeFraction
    rve4 = set_param(rve, p_e, 0.15)
    @test get_param(rve4, p_e) ≈ 0.15
    @test rve4.amounts[:C] isa CrackDensity

    # Original is untouched (no mutation)
    @test get_param(rve, p_f) ≈ 0.2

    # Setting matrix amount must error
    @test_throws ArgumentError set_param(rve, p_M, 0.5)
end

@testset "PropertyParameter — TensISO selectors" begin
    rve = RVE(:M)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:C => TensISO{3}(30.0, 10.0)))
    add_phase!(rve, :I, Ellipsoid(1.0), Dict(:C => TensISO{3}(60.0, 20.0));
               fraction = 0.2)

    # bulk / shear named selectors map to indices 1, 2 of TensISO{4,3}
    @test get_param(rve, property(:I, :C, :bulk))  ≈ 60.0
    @test get_param(rve, property(:I, :C, :shear)) ≈ 20.0
    @test get_param(rve, property(:I, :C, :K))     ≈ 60.0
    @test get_param(rve, property(:I, :C, :μ))     ≈ 20.0
    @test get_param(rve, property(:I, :C, 1))      ≈ 60.0
    @test get_param(rve, property(:I, :C, 2))      ≈ 20.0

    # set_param: bulk replaced, shear preserved
    rve2 = set_param(rve, property(:I, :C, :bulk), 90.0)
    @test get_param(rve2, property(:I, :C, :bulk))  ≈ 90.0
    @test get_param(rve2, property(:I, :C, :shear)) ≈ 20.0    # unchanged

    # Original tensor reference untouched
    @test get_param(rve, property(:I, :C, :bulk))  ≈ 60.0
end

@testset "PropertyParameter — error on unknown selector" begin
    rve = RVE(:M)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:C => TensISO{3}(30.0, 10.0)))
    add_phase!(rve, :I, Ellipsoid(1.0), Dict(:C => TensISO{3}(60.0, 20.0));
               fraction = 0.2)

    @test_throws ArgumentError get_param(rve, property(:I, :C, :nonexistent))
    @test_throws ArgumentError get_param(rve, property(:I, :K, :scalar))   # no :K key
end

@testset "GeometryParameter — Ellipsoid semi_axes" begin
    rve = RVE(:M)
    add_matrix!(rve, Ellipsoid(2.0, 1.5, 1.0),
                Dict(:C => TensISO{3}(30.0, 10.0)))
    add_phase!(rve, :I, Ellipsoid(1.0, 1.0, 0.5),
                Dict(:C => TensISO{3}(60.0, 20.0)); fraction = 0.2)

    p_a3 = geometry(:I, :semi_axes, 3)
    p_a1 = geometry(:I, :semi_axes, 1)

    @test get_param(rve, p_a3) ≈ 0.5
    @test get_param(rve, p_a1) ≈ 1.0

    # Update one axis, others preserved
    rve2 = set_param(rve, p_a3, 0.7)
    @test get_param(rve2, p_a3) ≈ 0.7
    @test get_param(rve2, p_a1) ≈ 1.0
    # Original untouched
    @test get_param(rve, p_a3) ≈ 0.5

    # Type-promotion to Dual via geometry param
    dx = ForwardDiff.Dual{Nothing, Float64, 1}(0.5, ForwardDiff.Partials((1.0,)))
    rve3 = set_param(rve, p_a3, dx)
    @test rve3.phases[:I].geometry.semi_axes[1] isa ForwardDiff.Dual
end

@testset "DistributionShapeParameter — UniformDistribution(Ellipsoid)" begin
    shape = Ellipsoid(2.0, 1.0, 0.5)
    rve = RVE(:M; distribution_shape = shape)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:C => TensISO{3}(30.0, 10.0)))
    add_phase!(rve, :I, Ellipsoid(1.0), Dict(:C => TensISO{3}(60.0, 20.0));
               fraction = 0.2)

    p_d2 = shape_param(:semi_axes, 2)
    @test get_param(rve, p_d2) ≈ 1.0

    rve2 = set_param(rve, p_d2, 1.7)
    @test get_param(rve2, p_d2) ≈ 1.7
    # Other axes preserved
    @test get_param(rve2, shape_param(:semi_axes, 1)) ≈ 2.0
    @test get_param(rve2, shape_param(:semi_axes, 3)) ≈ 0.5
end

@testset "_set_many — batch composition" begin
    rve = RVE(:M)
    add_matrix!(rve, Ellipsoid(1.0), Dict(:C => TensISO{3}(30.0, 10.0)))
    add_phase!(rve, :I, Ellipsoid(1.0), Dict(:C => TensISO{3}(60.0, 20.0));
               fraction = 0.2)

    ps = [amount(:I), property(:I, :C, :bulk), property(:M, :C, :shear)]
    vs = [0.3, 90.0, 12.0]

    # Internal helper used by gradient/jacobian
    rve2 = MeanFieldHom.Schemes._set_many(rve, ps, vs)
    @test get_param(rve2, ps[1]) ≈ 0.3
    @test get_param(rve2, ps[2]) ≈ 90.0
    @test get_param(rve2, ps[3]) ≈ 12.0
    # Original untouched
    @test get_param(rve, ps[1]) ≈ 0.2
end
