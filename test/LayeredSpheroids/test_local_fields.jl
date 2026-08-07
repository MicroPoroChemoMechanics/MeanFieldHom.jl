using Test
using MeanFieldHomogenization
using MeanFieldHomogenization.LayeredSpheroids: _base_fond
using TensND

# =============================================================================
#  test_local_fields.jl — pointwise temperature/gradient/flux
#  reconstruction (`localfields.jl`), checked against the two physical
#  continuity conditions at a perfect interface (normal flux and
#  tangential gradient continuous, NOT the full field vector) and the
#  far-field limit (recovers the imposed remote uniform gradient).
# =============================================================================

const K1 = TensISO{3}(5.0)
const K2 = TensISO{3}(20.0)
const K0 = TensISO{3}(2.0)

function _two_layer_spheroid(; Nseries = 8)
    focal2 = 8.0
    a_in = 2.9
    b_in = sqrt(a_in^2 - focal2)
    return LayeredSpheroid((a_in, 3.0), (b_in, 1.0), (K1, K2); Nseries)
end

@testset "local fields — continuity across a perfect interface" begin
    s = _two_layer_spheroid()
    q1 = real(s.q[1])
    p, φ = 0.4, 0.7
    ε = 1.0e-7

    e_q, e_p, e_φ = _base_fond(q1, p, φ)
    e_q, e_p, e_φ = real.(e_q), real.(e_p), real.(e_φ)

    u_below = local_flux(s, K0, q1 - ε, p, φ; H_axial = 1.0, H_trans = 0.3)
    u_above = local_flux(s, K0, q1 + ε, p, φ; H_axial = 1.0, H_trans = 0.3)
    un_below = sum(u_below .* e_q)
    un_above = sum(u_above .* e_q)
    @test un_below ≈ un_above rtol = 1.0e-4

    g_below = local_gradient(s, K0, q1 - ε, p, φ; H_axial = 1.0, H_trans = 0.3)
    g_above = local_gradient(s, K0, q1 + ε, p, φ; H_axial = 1.0, H_trans = 0.3)
    @test sum(g_below .* e_p) ≈ sum(g_above .* e_p) rtol = 1.0e-4
    @test sum(g_below .* e_φ) ≈ sum(g_above .* e_φ) rtol = 1.0e-4

    T_below = local_temperature(s, K0, q1 - ε, p, φ; H_axial = 1.0, H_trans = 0.3)
    T_above = local_temperature(s, K0, q1 + ε, p, φ; H_axial = 1.0, H_trans = 0.3)
    @test T_below ≈ T_above rtol = 1.0e-4
end

@testset "local fields — far field recovers the remote uniform gradient" begin
    s = _two_layer_spheroid()
    qbig = 1.0e3
    p, φ = 0.4, 0.7

    T = local_temperature(s, K0, qbig, p, φ; H_axial = 1.0, H_trans = 0.0)
    z = real(s.c * qbig * p)   # remote T ~ H₃·z, z = c·q·p (eq:xLeg, P₁(q)=q, P₁(p)=p)
    @test T ≈ z rtol = 1.0e-8

    g_axial = local_gradient(s, K0, qbig, p, φ; H_axial = 1.0, H_trans = 0.0)
    @test collect(g_axial) ≈ [0.0, 0.0, 1.0] atol = 1.0e-8

    g_trans = local_gradient(s, K0, qbig, 0.0, 0.0; H_axial = 0.0, H_trans = 1.0)
    @test collect(g_trans) ≈ [1.0, 0.0, 0.0] atol = 1.0e-8
end

@testset "local fields — accept both scalar and TensISO matrix conductivity" begin
    s = _two_layer_spheroid()
    T_tens = local_temperature(s, K0, 5.0, 0.2, 0.1)
    T_scalar = local_temperature(s, 2.0, 5.0, 0.2, 0.1)
    @test T_tens ≈ T_scalar
end
