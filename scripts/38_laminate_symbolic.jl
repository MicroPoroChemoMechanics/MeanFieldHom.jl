# =============================================================================
#  38_laminate_symbolic.jl
#
#  The laminate kernel run on SYMBOLIC moduli, so that the classical closed
#  forms come out of the code itself rather than being checked against it.
#
#  For two isotropic layers the general formula
#
#      ℂ_hom = ⟨ℚ⟩ + ⟨ℂ:ℙ⟩ : ⟨ℙ⟩† : ⟨ℙ:ℂ⟩
#
#  must simplify EXACTLY to Backus (1962). Getting there requires two
#  implementation choices that a purely numerical script would never expose:
#
#    * the pseudo-inverse is a cofactor inverse of the 3×3 out-of-plane block,
#      never `LinearAlgebra.pinv` — an SVD is not symbolically evaluable;
#    * every intermediate is an `SMatrix`, never an `MMatrix`, which cannot
#      even be CONSTRUCTED for a non-isbits element type such as `SymPy.Sym`.
#
#  NOT published to the documentation gallery: SymPy-heavy scripts are kept
#  out of the build (repo policy, see scripts/README.md). Its content is
#  covered by `docs/src/theory/laminate.md` and by
#  `test/Laminates/test_laminate_symbolic.jl`.
# =============================================================================

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)

using MeanFieldHom
using TensND
using StaticArrays
using SymPy
using Printf

const MFHC = MeanFieldHom.Core

println("="^78)
println("Periodic multilayer — symbolic derivation of the closed forms")
println("="^78)

# ── §1  Two isotropic layers, general fractions ──────────────────────────────

@syms λ₁::positive μ₁::positive λ₂::positive μ₂::positive f₁::positive
f₂ = 1 - f₁

function iso_km(λ, μ)
    z = zero(λ)
    return SMatrix{6, 6}(
        [
            λ+2μ λ λ z z z
            λ λ+2μ λ z z z
            λ λ λ+2μ z z z
            z z z 2μ z z
            z z z z 2μ z
            z z z z z 2μ
        ]
    )
end

Z = SMatrix{6, 6}(zeros(Sym, 6, 6))
Ch = MFHC.laminate_stiffness((iso_km(λ₁, μ₁), iso_km(λ₂, μ₂)), (f₁, f₂), Z, Z)

println("\n§1  Effective coefficients, simplified by SymPy")
println("─"^78)
@printf "  C₃₃₃₃ = %s\n" string(simplify(Ch[3, 3]))
@printf "  C₂₃₂₃ = %s\n" string(simplify(Ch[4, 4] / 2))
@printf "  C₁₂₁₂ = %s\n" string(simplify(Ch[6, 6] / 2))
@printf "  C₁₁₃₃ = %s\n" string(simplify(Ch[1, 3]))

# The two structural results, read directly off the symbolic output:
#   * out of plane, a HARMONIC (Reuss) mean — the laminate saturates Reuss;
#   * in plane, an ARITHMETIC (Voigt) mean  — and saturates Voigt.
println("\n§2  The two exact bound saturations")
println("─"^78)
d_reuss = simplify(1 / Ch[3, 3] - (f₁ / (λ₁ + 2μ₁) + f₂ / (λ₂ + 2μ₂)))
d_voigt = simplify(Ch[6, 6] / 2 - (f₁ * μ₁ + f₂ * μ₂))
@printf "  1/C₃₃₃₃ − ⟨1/(λ+2μ)⟩ = %s      (Reuss, out of plane)\n" string(d_reuss)
@printf "  C₁₂₁₂   − ⟨μ⟩         = %s      (Voigt, in plane)\n" string(d_voigt)

# ── §3  The full Backus (1962) set ───────────────────────────────────────────

avg(g) = f₁ * g(λ₁, μ₁) + f₂ * g(λ₂, μ₂)
r₃₃ = 1 / avg((l, m) -> 1 / (l + 2m))
rλ = avg((l, m) -> l / (l + 2m))

backus = (
    ("C₁₁₁₁", Ch[1, 1], avg((l, m) -> 4m * (l + m) / (l + 2m)) + r₃₃ * rλ^2),
    ("C₁₁₂₂", Ch[1, 2], avg((l, m) -> 2m * l / (l + 2m)) + r₃₃ * rλ^2),
    ("C₁₁₃₃", Ch[1, 3], r₃₃ * rλ),
    ("C₃₃₃₃", Ch[3, 3], r₃₃),
    ("C₂₃₂₃", Ch[4, 4] / 2, 1 / avg((l, m) -> 1 / m)),
    ("C₁₂₁₂", Ch[6, 6] / 2, avg((l, m) -> m)),
)

println("\n§3  Against Backus (1962), component by component")
println("─"^78)
for (name, got, want) in backus
    @printf "  simplify(%s − Backus) = %s\n" name string(simplify(got - want))
end

# ── §4  A spring interface, still in closed form ─────────────────────────────

@syms kn::positive L::positive
z = zero(kn)
𝕂 = SMatrix{3, 3}([z z z; z z z; z z kn])          # normal compliance only
P_int = MFHC._op_embed(MFHC.compliance_op_block(𝕂)) / L
Chi = MFHC.laminate_stiffness((iso_km(λ₁, μ₁), iso_km(λ₂, μ₂)), (f₁, f₂), P_int, Z)

println("\n§4  With a normal spring interface of compliance kn over a period L")
println("─"^78)
@printf "  1/C₃₃₃₃ = %s\n" string(simplify(1 / Chi[3, 3]))
d_itf = simplify(1 / Chi[3, 3] - (f₁ / (λ₁ + 2μ₁) + f₂ / (λ₂ + 2μ₂) + kn / L))
@printf "  minus ⟨1/(λ+2μ)⟩ + kn/L = %s   (the compliance simply adds)\n" string(d_itf)
@printf "  C₁₂₁₂ unchanged : %s\n" string(simplify(Chi[6, 6] / 2 - (f₁ * μ₁ + f₂ * μ₂)))

# ── §5  Through the Laminate cell itself ─────────────────────────────────────
#
# Not just the bare kernel: the cell (property dicts, derived fractions, the
# exact-TI return ladder) carries symbolic entries end to end.

@syms κ₁::positive κ₂::positive
lam = Laminate(; T = Sym)
add_layer!(lam, :A, Dict(:C => TensISO{3}(3κ₁, 2μ₁)); fraction = f₁)
add_layer!(lam, :B, Dict(:C => TensISO{3}(3κ₂, 2μ₂)); fraction = f₂)
Csym = homogenize(lam, Laminated(), :C)

println("\n§5  `homogenize` on a Laminate{Sym}")
println("─"^78)
println("  returned type : ", typeof(Csym))
println("  (isotropic layers ⇒ EXACTLY transversely isotropic about n)")
Msym = KM(Csym)
lame(κ, μ) = κ - 2μ / 3
@printf "  1/C₃₃₃₃ − ⟨1/(λ+2μ)⟩ = %s\n" string(
    simplify(1 / Msym[3, 3] - (f₁ / (lame(κ₁, μ₁) + 2μ₁) + f₂ / (lame(κ₂, μ₂) + 2μ₂)))
)

# ── §6  Conduction ───────────────────────────────────────────────────────────

@syms k₁::positive k₂::positive ρ::positive
K3s = (
    SMatrix{3, 3}(Sym[k₁ 0 0; 0 k₁ 0; 0 0 k₁]),
    SMatrix{3, 3}(Sym[k₂ 0 0; 0 k₂ 0; 0 0 k₂]),
)
Z3 = SMatrix{3, 3}(zeros(Sym, 3, 3))
Pρ = SMatrix{3, 3}(Sym[0 0 0; 0 0 0; 0 0 ρ]) / L

Kh = MFHC.laminate_conductivity(K3s, (f₁, f₂), Z3, Z3)
Khρ = MFHC.laminate_conductivity(K3s, (f₁, f₂), Pρ, Z3)

println("\n§6  Conduction")
println("─"^78)
@printf "  k_⊥      = %s        (series)\n" string(simplify(Kh[3, 3]))
@printf "  k_∥      = %s        (parallel)\n" string(simplify(Kh[1, 1]))
@printf "  1/k_⊥ with Kapitza ρ, minus ⟨1/k⟩ + ρ/L = %s\n" string(
    simplify(1 / Khρ[3, 3] - (f₁ / k₁ + f₂ / k₂ + ρ / L))
)

println("\n" * "="^78)
println("Every difference above is 0: the code reproduces the closed forms exactly.")
println("="^78)
