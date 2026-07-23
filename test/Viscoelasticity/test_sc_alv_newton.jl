using Test
using MeanFieldHom
using TensND
using LinearAlgebra

# =============================================================================
#  test_sc_alv_newton.jl — Newton-Raphson ligne par ligne pour le SC ALV
#  (`src/Viscoelasticity/schemes_alv_sc_newton.jl`).
#
#  NOTE : `self_consistent_alv_newton` n'est ni exporté ni appelé depuis le
#  reste du paquet — il s'atteint uniquement par son chemin qualifié. On le
#  teste ici contre les deux références disponibles : le SC élastique dans la
#  limite de Heaviside, et le SC ALV Anderson-Picard en viscoélastique.
# =============================================================================

const _sc_newton = MeanFieldHom.Viscoelasticity.self_consistent_alv_newton
const _to_mandel66 = MeanFieldHom.Viscoelasticity._tens_to_mandel66

@testset "sc_alv_newton — limite élastique (Heaviside)" begin
    rve = RVE(:M)
    add_matrix!(
        rve, Ellipsoid(1.0, 1.0, 1.0),
        Dict(:C => heaviside_law(TensISO{3}(30.0, 8.0)))
    )
    add_phase!(
        rve, :I, Ellipsoid(1.0, 1.0, 1.0),
        Dict(:C => heaviside_law(TensISO{3}(60.0, 16.0)));
        fraction = 0.2
    )

    times = collect(0.0:0.5:1.5)
    n = length(times)
    C_newton = _sc_newton(rve, :C; times = times, abstol = 1.0e-12)

    @test size(C_newton) == (6n, 6n)

    # Référence : SC élastique sur les mêmes modules.
    rve_e = RVE(:M)
    add_matrix!(rve_e, Ellipsoid(1.0, 1.0, 1.0), Dict(:C => TensISO{3}(30.0, 8.0)))
    add_phase!(
        rve_e, :I, Ellipsoid(1.0, 1.0, 1.0),
        Dict(:C => TensISO{3}(60.0, 16.0)); fraction = 0.2
    )
    C_e = _to_mandel66(homogenize(rve_e, SelfConsistent(), :C))

    # En Heaviside la réponse est constante : blocs diagonaux = SC élastique,
    # blocs hors diagonale nuls.  Le Newton ligne par ligne et le SC élastique
    # sont deux solveurs distincts : ils convergent vers la même racine à leur
    # tolérance propre près (~1e-10 en relatif), d'où un `rtol` plutôt qu'un
    # `atol` serré.
    for i in 1:n
        rows = (6 * (i - 1) + 1):(6 * i)
        @test isapprox(C_newton[rows, rows], C_e; rtol = 1.0e-7, atol = 1.0e-10)
        for j in 1:(i - 1)
            cols = (6 * (j - 1) + 1):(6 * j)
            @test maximum(abs, C_newton[rows, cols]) ≤ 1.0e-10
        end
    end
end

@testset "sc_alv_newton — accord avec le SC ALV Anderson-Picard" begin
    # Les deux solveurs cherchent la même racine `C_eff = step(C_eff)` ; sur
    # une configuration à faible contraste ils doivent tomber sur le même
    # point fixe.
    rve = RVE(:M)
    add_matrix!(
        rve, Ellipsoid(1.0, 1.0, 1.0),
        Dict(:C => maxwell_iso(10.0, 4.0, 1.0, 0.5))
    )
    add_phase!(
        rve, :I, Ellipsoid(1.0, 1.0, 1.0),
        Dict(:C => heaviside_law(TensISO{3}(20.0, 8.0)));
        fraction = 0.15
    )

    times = collect(0.0:0.25:1.0)
    C_newton = _sc_newton(rve, :C; times = times, abstol = 1.0e-12)
    C_picard = self_consistent_alv(
        rve, :C; times = times, abstol = 1.0e-12,
        maxiters = 2000
    )

    @test isapprox(C_newton, C_picard; rtol = 1.0e-6, atol = 1.0e-9)
end

@testset "sc_alv_newton — causalité (triangulaire inférieure par blocs)" begin
    rve = RVE(:M)
    add_matrix!(
        rve, Ellipsoid(1.0, 1.0, 1.0),
        Dict(:C => maxwell_iso(12.0, 5.0, 1.0, 0.4))
    )
    add_phase!(
        rve, :I, Ellipsoid(1.0, 1.0, 1.0),
        Dict(:C => heaviside_law(TensISO{3}(30.0, 12.0)));
        fraction = 0.1
    )

    times = collect(0.0:0.25:0.75)
    n = length(times)
    C = _sc_newton(rve, :C; times = times, abstol = 1.0e-12)

    # L'équation SC ALV est causale : la ligne i ne dépend que des lignes ≤ i,
    # donc le bloc supérieur strict doit être nul.
    for i in 1:n, j in (i + 1):n
        rows = (6 * (i - 1) + 1):(6 * i)
        cols = (6 * (j - 1) + 1):(6 * j)
        @test maximum(abs, C[rows, cols]) ≤ 1.0e-10
    end
end

@testset "sc_alv_newton — phase fissurée (CrackDensity)" begin
    rve = RVE(:M)
    add_matrix!(
        rve, Ellipsoid(1.0, 1.0, 1.0),
        Dict(:C => maxwell_iso(10.0, 4.0, 1.0, 0.5))
    )
    add_phase!(
        rve, :F, Ellipsoid(1.0, 1.0, 0.0),
        Dict(:C => heaviside_law(TensISO{3}(1.0e-9, 1.0e-9)));
        density = 0.05
    )

    times = collect(0.0:0.25:0.75)
    n = length(times)
    C = _sc_newton(rve, :C; times = times, abstol = 1.0e-10)

    @test size(C) == (6n, 6n)
    @test all(isfinite, C)

    # Les fissures assouplissent : à t = 0 la rigidité effective doit être
    # strictement inférieure à celle de la matrice seule.
    C_M = MeanFieldHom.Viscoelasticity._trapezoidal_relaxation(
        matrix_property(rve, :C), times, 6
    )
    @test C[1, 1] < C_M[1, 1]
end

@testset "sc_alv_newton — erreurs d'entrée" begin
    # Propriété de matrice élastique (pas une ViscoLaw).
    rve_e = RVE(:M)
    add_matrix!(rve_e, Ellipsoid(1.0, 1.0, 1.0), Dict(:C => TensISO{3}(30.0, 8.0)))
    @test_throws ArgumentError _sc_newton(rve_e, :C; times = [0.0, 1.0])

    # Propriété de phase élastique alors que la matrice est visqueuse.
    rve_p = RVE(:M)
    add_matrix!(
        rve_p, Ellipsoid(1.0, 1.0, 1.0),
        Dict(:C => heaviside_law(TensISO{3}(30.0, 8.0)))
    )
    add_phase!(
        rve_p, :I, Ellipsoid(1.0, 1.0, 1.0),
        Dict(:C => TensISO{3}(60.0, 16.0)); fraction = 0.1
    )
    @test_throws ArgumentError _sc_newton(rve_p, :C; times = [0.0, 1.0])

    # CrackDensity sur une géométrie qui n'est pas une fissure.
    rve_c = RVE(:M)
    add_matrix!(
        rve_c, Ellipsoid(1.0, 1.0, 1.0),
        Dict(:C => heaviside_law(TensISO{3}(30.0, 8.0)))
    )
    add_phase!(
        rve_c, :F, Ellipsoid(1.0, 1.0, 1.0),
        Dict(:C => heaviside_law(TensISO{3}(1.0e-9, 1.0e-9)));
        density = 0.05
    )
    @test_throws ArgumentError _sc_newton(rve_c, :C; times = [0.0, 1.0])
end
