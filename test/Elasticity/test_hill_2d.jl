using Test
using MeanFieldHom
using TensND
using LinearAlgebra
import MeanFieldHom.Elasticity: _hill_2d_iso, _hill_2d_aniso

# =============================================================================
#  test_hill_2d.jl — tenseur de Hill 2D, forme fermée isotrope et quadrature
#  anisotrope.
#
#  Les deux chemins sont atteints depuis le MÊME appel `hill_tensor(ell, C₀)` :
#  `_hill_2d_iso` quand `C₀` est un `TensISO{4,2}`, `_hill_2d_aniso` sinon.
#  Ils doivent donc coïncider quand `C₀` est isotrope mais stocké sous forme
#  générique — c'est le verrou principal de ce fichier.
#
#  Référence externe : Mura (1987) eq. 11.22, tenseur d'Eshelby du cylindre
#  elliptique, avec P = S : C₀⁻¹.  Convention du paquet :
#  C₀ = TensISO{2}(α, β), α = 3k, β = 2μ, déformation plane ⇒ ν = (α-β)/(2α).
# =============================================================================

# Même tenseur, exprimé en `Tens` générique : force le routage anisotrope.
function _as_generic(C)
    A = zeros(2, 2, 2, 2)
    for i in 1:2, j in 1:2, p in 1:2, q in 1:2
        A[i, j, p, q] = C[i, j, p, q]
    end
    return TensND.Tens(A)
end

# Mura (1987) eq. 11.22 — cylindre elliptique, demi-axes (a, b), déformation
# plane, coefficient de Poisson ν.
function _mura_S(a, b, nu)
    s = 1 / (2 * (1 - nu))
    ab = (a + b)^2
    S = zeros(2, 2, 2, 2)
    S[1, 1, 1, 1] = s * ((b^2 + 2a * b) / ab + (1 - 2nu) * b / (a + b))
    S[2, 2, 2, 2] = s * ((a^2 + 2a * b) / ab + (1 - 2nu) * a / (a + b))
    S[1, 1, 2, 2] = s * (b^2 / ab - (1 - 2nu) * b / (a + b))
    S[2, 2, 1, 1] = s * (a^2 / ab - (1 - 2nu) * a / (a + b))
    v = s * ((a^2 + b^2) / (2ab) + (1 - 2nu) / 2)
    for (i, j, p, q) in ((1, 2, 1, 2), (1, 2, 2, 1), (2, 1, 1, 2), (2, 1, 2, 1))
        S[i, j, p, q] = v
    end
    return S
end

function _mura_P(a, b, k, mu)
    α, β = 3k, 2mu
    nu = (α - β) / (2α)
    S = _mura_S(a, b, nu)
    Cinv = TensISO{2}(1 / α, 1 / β)
    return [
        sum(S[i, j, m, n] * Cinv[m, n, p, q] for m in 1:2, n in 1:2)
            for i in 1:2, j in 1:2, p in 1:2, q in 1:2
    ]
end

_maxdiff(P, Q) = maximum(abs(P[i, j, p, q] - Q[i, j, p, q]) for i in 1:2, j in 1:2, p in 1:2, q in 1:2)

@testset "hill 2D — forme fermée isotrope vs Mura (1987)" begin
    for (k, mu) in ((5.0, 2.0), (1.0, 1.0), (10.0, 0.5), (0.5, 4.0))
        for rho in (1.0, 0.8, 0.5, 0.2)
            ell = Ellipsoid(1.0, rho)
            P = hill_tensor(ell, TensISO{2}(3k, 2mu))
            @test _maxdiff(P, _mura_P(1.0, rho, k, mu)) < 1.0e-12
        end
    end
end

@testset "hill 2D — quadrature anisotrope vs Mura (1987)" begin
    for (k, mu) in ((5.0, 2.0), (1.0, 1.0)), rho in (1.0, 0.5, 0.3)
        ell = Ellipsoid(1.0, rho)
        P = _hill_2d_aniso(ell, _as_generic(TensISO{2}(3k, 2mu)))
        @test _maxdiff(P, _mura_P(1.0, rho, k, mu)) < 1.0e-9
    end
end

@testset "hill 2D — les deux chemins coïncident sur le même C₀" begin
    # Régression : la forme fermée `_hill_2d_iso` et la quadrature générale
    # `_hill_2d_aniso` partent du même appel `hill_tensor` et doivent donner
    # le même tenseur.  Elles divergeaient de ~1e-2 avant correction de la
    # forme fermée.
    for (k, mu) in ((5.0, 2.0), (2.0, 3.0), (10.0, 0.5)), rho in (1.0, 0.7, 0.4)
        ell = Ellipsoid(1.0, rho)
        C = TensISO{2}(3k, 2mu)
        @test _maxdiff(hill_tensor(ell, C), _hill_2d_aniso(ell, _as_generic(C))) < 1.0e-8
    end
end

@testset "hill 2D — symétries et signe" begin
    for rho in (1.0, 0.6)
        P = hill_tensor(Ellipsoid(1.0, rho), TensISO{2}(15.0, 4.0))
        # Symétrie majeure et symétries mineures.
        for i in 1:2, j in 1:2, p in 1:2, q in 1:2
            @test P[i, j, p, q] ≈ P[p, q, i, j] atol = 1.0e-12
            @test P[i, j, p, q] ≈ P[j, i, p, q] atol = 1.0e-12
            @test P[i, j, p, q] ≈ P[i, j, q, p] atol = 1.0e-12
        end
        @test P[1, 1, 1, 1] > 0
        @test P[2, 2, 2, 2] > 0
        @test P[1, 2, 1, 2] > 0
    end
end

@testset "hill 2D — limite incompressible (k = Inf)" begin
    mu = 2.0
    for rho in (1.0, 0.5, 0.25)
        ell = Ellipsoid(1.0, rho)
        P_inf = hill_tensor(ell, TensISO{2}(Inf, 2mu))
        P_big = hill_tensor(ell, TensISO{2}(3.0e12, 2mu))
        @test all(
            isfinite(P_inf[i, j, p, q])
                for i in 1:2, j in 1:2, p in 1:2, q in 1:2
        )
        # La branche `isinf` doit être la limite continue de la branche
        # générale.
        @test _maxdiff(P_inf, P_big) < 1.0e-9
        # Incompressibilité : la partie sphérique de P s'annule.
        @test abs(P_inf[1, 1, 1, 1] + P_inf[1, 1, 2, 2]) < 1.0e-12
    end
end

@testset "hill 2D — cercle : P est isotrope 2D" begin
    k, mu = 5.0, 2.0
    P = hill_tensor(Ellipsoid(1.0, 1.0), TensISO{2}(3k, 2mu))
    @test P isa TensND.TensISO{4, 2}
    # Valeurs propres fermées : P_J = 1/(3k+2μ), P_K = (3k+4μ)/(4μ(3k+2μ)).
    PJ, PK = TensND.get_data(P)
    @test PJ ≈ 1 / (3k + 2mu)
    @test PK ≈ (3k + 4mu) / (4mu * (3k + 2mu))
end
