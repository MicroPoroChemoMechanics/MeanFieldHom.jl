using Test
using MeanFieldHom
using TensND
using LinearAlgebra
import MeanFieldHom.Viscoelasticity: _is_ortho_block, _is_ti_block, _is_iso_block,
    _try_iso_pairs, _try_ti_tuples, _try_ortho_tuples, _ortho_pair, _ortho_blocks

# =============================================================================
#  test_ortho_dispatch_alv.jl — chemins rapides ORTHO de `_homogenize_alv_dispatch`
#  (`src/Viscoelasticity/homogenize_alv.jl`).
#
#  `test_ortho_alv.jl` teste les primitives ortho (`voigt_alv_ortho`, …) en les
#  appelant directement, mais tous ses RVE sont isotropes : comme iso ⊂ TI ⊂
#  ortho, le dispatcher prend alors le raccourci iso et les branches ortho de
#  `_homogenize_alv_dispatch` ne sont jamais exécutées.
#
#  Ici on construit un RVE franchement orthotrope (ni iso ni TI) pour forcer
#  `_try_iso_pairs` et `_try_ti_tuples` à rendre `nothing` et faire tomber le
#  dispatcher sur `_try_ortho_tuples`.
# =============================================================================

const _ORTHO_FRAME = TensND.CanonicalBasis{3, Float64}()

# Tenseur orthotrope franc : les trois modules normaux et les trois modules de
# cisaillement sont tous distincts, donc ni isotrope ni TI autour d'aucun axe.
_ortho_tensor(s) = TensND.TensOrtho(
    20.0s, 8.0s, 6.0s,       # C11, C12, C13
    30.0s, 7.0s, 40.0s,      # C22, C23, C33
    5.0s, 6.0s, 7.0s,        # C44, C55, C66
    _ORTHO_FRAME
)

function _ortho_rve(; fraction = 0.25)
    rve = RVE(:M)
    add_matrix!(
        rve, Ellipsoid(1.0, 1.0, 1.0),
        Dict(:C => heaviside_law(_ortho_tensor(1.0)))
    )
    add_phase!(
        rve, :I, Ellipsoid(1.0, 1.0, 1.0),
        Dict(:C => heaviside_law(_ortho_tensor(2.0)));
        fraction = fraction
    )
    return rve
end

@testset "ortho dispatch — le RVE est bien hors des raccourcis iso et TI" begin
    times = collect(0.0:0.5:1.0)
    M_M = trapezoidal_matrix(heaviside_law(_ortho_tensor(1.0)), times)
    M_I = trapezoidal_matrix(heaviside_law(_ortho_tensor(2.0)), times)

    # Prérequis du test : sans ça le dispatcher partirait sur le chemin iso ou
    # TI et les branches visées resteraient mortes.
    @test !_is_iso_block(M_M)
    @test !_is_ti_block(M_M)
    @test _is_ortho_block(M_M)

    @test _try_iso_pairs([M_M, M_I]) === nothing
    @test _try_ti_tuples([M_M, M_I]) === nothing

    o = _try_ortho_tuples([M_M, M_I])
    @test o !== nothing
    @test length(o) == 2
    @test length(o[1]) == 12

    # `_try_ortho_tuples` doit rendre `nothing` dès qu'une seule matrice sort
    # de la forme ortho.
    M_bad = copy(M_M)
    M_bad[1, 4] += 1.0
    @test !_is_ortho_block(M_bad)
    @test _try_ortho_tuples([M_M, M_bad]) === nothing

    # Cas vide : vecteur vide, pas `nothing`.
    empty_in = Matrix{Float64}[]
    @test _try_ortho_tuples(empty_in) == NTuple{12, Matrix{Float64}}[]
    @test _try_ti_tuples(empty_in) == NTuple{6, Matrix{Float64}}[]
    @test _try_iso_pairs(empty_in) == Tuple{Matrix{Float64}, Matrix{Float64}}[]
end

@testset "ortho dispatch — Voigt et Reuss passent par le chemin ortho" begin
    times = collect(0.0:0.5:1.5)
    n = length(times)
    f = 0.25
    rve = _ortho_rve(fraction = f)

    M_M = trapezoidal_matrix(heaviside_law(_ortho_tensor(1.0)), times)
    M_I = trapezoidal_matrix(heaviside_law(_ortho_tensor(2.0)), times)
    fr = [1 - f, f]

    C_voigt = homogenize_alv(rve, Voigt(), :C; times = times)
    @test size(C_voigt) == (6n, 6n)
    @test _is_ortho_block(C_voigt)
    @test !_is_ti_block(C_voigt)                 # on est bien resté ortho
    @test isapprox(C_voigt, voigt_alv([M_M, M_I], fr); atol = 1.0e-10)

    C_reuss = homogenize_alv(rve, Reuss(), :C; times = times)
    @test _is_ortho_block(C_reuss)
    @test isapprox(C_reuss, reuss_alv([M_M, M_I], fr); atol = 1.0e-8)

    # Encadrement de Voigt-Reuss sur les termes diagonaux.
    for i in 1:(6n)
        @test C_reuss[i, i] ≤ C_voigt[i, i] + 1.0e-9
    end
end

@testset "ortho dispatch — Dilute, Mori-Tanaka et Maxwell restent ortho" begin
    times = collect(0.0:0.5:1.5)
    n = length(times)
    rve = _ortho_rve(fraction = 0.2)

    for scheme in (Dilute(), MoriTanaka(), Maxwell())
        C = homogenize_alv(rve, scheme, :C; times = times)
        @test size(C) == (6n, 6n)
        @test all(isfinite, C)
        # La forme ortho est stable par ces schémas : phases ortho coaxiales
        # ⇒ résultat ortho dans le même repère matériel.
        @test _is_ortho_block(C)
        @test !_is_iso_block(C)
    end
end

@testset "ortho dispatch — fraction nulle redonne la matrice" begin
    times = collect(0.0:0.5:1.0)
    rve = _ortho_rve(fraction = 0.0)
    M_M = trapezoidal_matrix(heaviside_law(_ortho_tensor(1.0)), times)

    for scheme in (Voigt(), Reuss(), Dilute(), MoriTanaka())
        C = homogenize_alv(rve, scheme, :C; times = times)
        @test isapprox(C, M_M; atol = 1.0e-8)
    end
end

@testset "ortho dispatch — round-trip _ortho_pair / _ortho_blocks" begin
    times = collect(0.0:0.5:1.5)
    M = trapezoidal_matrix(heaviside_law(_ortho_tensor(1.3)), times)

    o = _ortho_pair(M)
    @test length(o) == 12
    @test isapprox(_ortho_blocks(o), M; atol = 1.0e-12)
end
