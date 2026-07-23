# =============================================================================
#  test_param_conversions.jl — Interprétations physiques des coefficients de
#  classe de symétrie (`src/Elasticity/param_conversions.jl`).
#
#  Couvre les deux familles exportées :
#    • isotrope   : k_mu / iso_stiffness, E_nu / iso_stiffness_E_nu
#    • TI (Hoenig): hoenig_params / hoenig_stiffness, y compris la forme à
#      deux arguments qui projette d'abord sur le span TI.
#
#  Chaque paire est vérifiée dans les deux sens (aller-retour), contre les
#  formules fermées de l'élasticité isotrope, et sur le cas dégénéré d'un
#  tenseur TI qui est en fait isotrope.
# =============================================================================

using Test
using MeanFieldHom
using TensND
using LinearAlgebra

@testset "k_mu / iso_stiffness — aller-retour" begin
    k, mu = 30.0, 12.0
    C = iso_stiffness(k, mu)
    @test C isa TensND.TensISO{4}

    # C = 3k·𝕁 + 2μ·𝕂 : les coefficients bruts sont (3k, 2μ).
    @test collect(TensND.get_data(C)) ≈ [3k, 2mu]

    k_back, mu_back = k_mu(C)
    @test k_back ≈ k
    @test mu_back ≈ mu

    # La souplesse S = C⁻¹ se lit avec le même appel, via inv (cf. l'en-tête
    # de param_conversions.jl : pas de variante « _compliance »).
    S = inv(C)
    k_S, mu_S = k_mu(inv(S))
    @test k_S ≈ k
    @test mu_S ≈ mu
end

@testset "E_nu / iso_stiffness_E_nu — aller-retour et cohérence avec (k, μ)" begin
    E, nu = 210.0, 0.3
    C = iso_stiffness_E_nu(E, nu)
    @test C isa TensND.TensISO{4}

    E_back, nu_back = E_nu(C)
    @test E_back ≈ E
    @test nu_back ≈ nu

    # Cohérence avec les relations fermées k = E/(3(1-2ν)), μ = E/(2(1+ν)).
    k, mu = k_mu(C)
    @test k ≈ E / (3 * (1 - 2nu))
    @test mu ≈ E / (2 * (1 + nu))

    # Et dans l'autre sens depuis (k, μ).
    E2, nu2 = E_nu(iso_stiffness(k, mu))
    @test E2 ≈ E
    @test nu2 ≈ nu

    # ν = 0 : k = E/3 et μ = E/2.
    C0 = iso_stiffness_E_nu(100.0, 0.0)
    k0, mu0 = k_mu(C0)
    @test k0 ≈ 100.0 / 3
    @test mu0 ≈ 50.0
    @test all(isapprox.(E_nu(C0), (100.0, 0.0); atol = 1.0e-12))
end

@testset "hoenig_params / hoenig_stiffness — aller-retour" begin
    axis = [0.0, 0.0, 1.0]
    E1, h, nu1, nu2, gamma = 12.0, 2.5, 0.25, 0.2, 1.4

    C = hoenig_stiffness(E1, h, nu1, nu2, gamma, axis)
    @test C isa TensND.TensTI{4}

    p = hoenig_params(C)
    @test p.E1 ≈ E1
    @test p.h ≈ h
    @test p.nu1 ≈ nu1
    @test p.nu2 ≈ nu2
    @test p.gamma ≈ gamma

    # Les paramètres sont nommés : l'ordre documenté doit être respecté.
    @test collect(p) ≈ [E1, h, nu1, nu2, gamma]

    # Aller-retour complet sur le tenseur lui-même, pas seulement sur les
    # paramètres.
    C_back = hoenig_stiffness(p.E1, p.h, p.nu1, p.nu2, p.gamma, axis)
    @test Array(C_back) ≈ Array(C)
end

@testset "hoenig_params — axe non canonique" begin
    axis = normalize([1.0, 1.0, 0.0])
    E1, h, nu1, nu2, gamma = 20.0, 1.8, 0.3, 0.15, 0.9

    C = hoenig_stiffness(E1, h, nu1, nu2, gamma, axis)
    p = hoenig_params(C)

    # Les paramètres de Hoenig sont intrinsèques : ils ne dépendent pas de
    # l'orientation de l'axe de symétrie.
    @test p.E1 ≈ E1
    @test p.h ≈ h
    @test p.nu1 ≈ nu1
    @test p.nu2 ≈ nu2
    @test p.gamma ≈ gamma
end

@testset "hoenig_params — forme à deux arguments (projection préalable)" begin
    axis = [0.0, 0.0, 1.0]
    E1, h, nu1, nu2, gamma = 15.0, 2.0, 0.28, 0.18, 1.1

    C_ti = hoenig_stiffness(E1, h, nu1, nu2, gamma, axis)

    # Un TensTI déjà dans le span TI doit se projeter sur lui-même : la forme
    # à deux arguments doit alors redonner exactement les mêmes paramètres
    # que la forme à un argument.
    C_full = Tens(Array(C_ti))
    p2 = hoenig_params(C_full, axis)
    p1 = hoenig_params(C_ti)

    @test p2.E1 ≈ p1.E1
    @test p2.h ≈ p1.h
    @test p2.nu1 ≈ p1.nu1
    @test p2.nu2 ≈ p1.nu2
    @test p2.gamma ≈ p1.gamma
end

@testset "hoenig_params — cas dégénéré isotrope" begin
    # Une TI construite depuis un tenseur isotrope doit rendre γ = 1 (pas
    # d'anisotropie de cisaillement) et retomber sur le ν isotrope.
    E, nu = 70.0, 0.25
    C_iso = iso_stiffness_E_nu(E, nu)
    axis = [0.0, 0.0, 1.0]

    C_ti, = TensND.proj_tens(Val(:TI), Array(C_iso), axis)
    p = hoenig_params(C_ti)

    @test p.gamma ≈ 1.0
    @test p.h ≈ 1.0
    @test p.nu1 ≈ nu
    @test p.nu2 ≈ nu
    @test p.E1 ≈ E
end

@testset "param_conversions — souplesse TI via inv" begin
    axis = [0.0, 0.0, 1.0]
    C = hoenig_stiffness(18.0, 2.2, 0.26, 0.17, 1.25, axis)
    S = inv(C)

    # inv sur un TensTI reste dans la classe TI, donc hoenig_params(inv(S))
    # doit redonner les paramètres de C.
    p = hoenig_params(inv(S))
    q = hoenig_params(C)
    @test p.E1 ≈ q.E1
    @test p.h ≈ q.h
    @test p.nu1 ≈ q.nu1
    @test p.nu2 ≈ q.nu2
    @test p.gamma ≈ q.gamma
end
