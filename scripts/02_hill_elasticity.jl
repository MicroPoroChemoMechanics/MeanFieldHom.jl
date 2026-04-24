# =============================================================================
#  02_hill_elasticity.jl
#
#  4th-order Hill polarization tensors for elastic inclusions.
#
#  Run from the MeanFieldHom.jl root:
#    julia --project=. scripts/02_hill_elasticity.jl
#
#  Sections:
#   § 1  Isotropic matrix — sphere, prolate, oblate, triaxial
#   § 2  Anisotropic matrix — residue vs DECUHR comparison
#   § 3  2D plane strain — circle, ellipse
#   § 4  Eshelby tensor  S^E = P : C₀
#   § 5  Dilute homogenization — porous solid
# =============================================================================

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io=devnull)

using MeanFieldHom
using TensND
using LinearAlgebra
using Printf

# ─── Reference material: isotropic steel  E = 210 GPa, ν = 0.3 ──────────────
const E_ref  = 210e3   # MPa
const ν_ref  = 0.3
const λ_ref  = E_ref * ν_ref / ((1 + ν_ref) * (1 - 2ν_ref))
const μ_ref  = E_ref / (2 * (1 + ν_ref))
const k_ref  = λ_ref + 2μ_ref / 3

const C_iso = TensISO{3}(3k_ref, 2μ_ref)

@printf "Reference matrix:  E = %.1f MPa,  ν = %.2f\n" E_ref ν_ref
@printf "  λ = %.2f,  μ = %.2f,  k = %.2f MPa\n\n" λ_ref μ_ref k_ref

# ─── Helpers ─────────────────────────────────────────────────────────────────

const voigt_idx = ((1,1),(2,2),(3,3),(2,3),(1,3),(1,2))
const voigt_lab = ["11","22","33","23","13","12"]

function print_voigt(C; label="P", scale=1e6, unit="×10⁻⁶ MPa⁻¹")
    println("  Voigt[$label] ($unit):")
    print("      ")
    for l in voigt_lab; @printf "%10s" l; end; println()
    for (I,(i,j)) in enumerate(voigt_idx)
        @printf "  %2s | " voigt_lab[I]
        for (k,l) in voigt_idx; @printf "%9.4f " scale*C[i,j,k,l]; end
        println()
    end
end

# Double contraction  (A:B)_{ijkl} = Σ_{mn} A_{ijmn} B_{mnkl}
function dcontract(A, B, dim=3)
    S = zeros(dim,dim,dim,dim)
    for i in 1:dim, j in 1:dim, k in 1:dim, l in 1:dim
        for m in 1:dim, n in 1:dim
            S[i,j,k,l] += A[i,j,m,n] * B[m,n,k,l]
        end
    end
    return S
end

# Mandel 6×6 matrix (M_{IJ} = T_{ijkl} fI fJ, f=[1,1,1,√2,√2,√2])
const mandel_f = (1.0,1.0,1.0,√2,√2,√2)
function mandel(C)
    [C[i,j,k,l] * mandel_f[I] * mandel_f[J]
     for (I,(i,j)) in enumerate(voigt_idx), (J,(k,l)) in enumerate(voigt_idx)]
end

# ═══════════════════════════════════════════════════════════════════════════
println("="^70)
println("  § 1  ISOTROPIC MATRIX — 3D geometries")
println("="^70)

# ─── 1a. Sphere ──────────────────────────────────────────────────────────────
println("\n── Sphere  (a = b = c = 1) ──────────────────────────────────────────")
let ell = Ellipsoid(1.0)
    P = hill_tensor(ell, C_iso)

    P1111_theory = 1/(5*(λ_ref+2μ_ref)) + (1/3 - 1/5)/μ_ref

    @printf "\n  P[1,1,1,1] = %12.9e MPa⁻¹\n" P[1,1,1,1]
    @printf "  Theory     = %12.9e MPa⁻¹  (err=%.2e)\n" P1111_theory abs(P[1,1,1,1]-P1111_theory)
    @printf "  Isotropy:  P[1111]=P[2222]? err=%.2e\n" abs(P[1,1,1,1]-P[2,2,2,2])

    println("\n  Voigt matrix:")
    print_voigt(P)
end

# ─── 1b. Prolate spheroid (fiber) ────────────────────────────────────────────
println("\n── Prolate spheroid  (a = 5, b = c = 1)  — long fiber ───────────────")
let ell = Ellipsoid(5.0, 1.0, 1.0)
    P = hill_tensor(ell, C_iso)

    @printf "\n  P[1,1,1,1] = %12.9e MPa⁻¹  (axial)\n" P[1,1,1,1]
    @printf "  P[2,2,2,2] = %12.9e MPa⁻¹  (transverse)\n" P[2,2,2,2]
    @printf "  P[3,3,3,3] = %12.9e MPa⁻¹\n" P[3,3,3,3]
    @printf "  Transverse isotropy: P[2222]=P[3333]? err=%.2e\n" abs(P[2,2,2,2]-P[3,3,3,3])
    println("\n  Voigt matrix:")
    print_voigt(P; label="P prolate a/b=5")
end

# ─── 1c. Oblate spheroid (disk) ──────────────────────────────────────────────
println("\n── Oblate spheroid  (a = b = 5, c = 1)  — disk / platelet ──────────")
let ell = Ellipsoid(5.0, 5.0, 1.0)
    P = hill_tensor(ell, C_iso)

    @printf "\n  P[1,1,1,1] = %12.9e MPa⁻¹  (in-plane)\n" P[1,1,1,1]
    @printf "  P[3,3,3,3] = %12.9e MPa⁻¹  (normal to disk — dominant)\n" P[3,3,3,3]
    @printf "  P[1111]=P[2222]? err=%.2e\n" abs(P[1,1,1,1]-P[2,2,2,2])
end

# ─── 1d. Triaxial ────────────────────────────────────────────────────────────
println("\n── Triaxial  (a = 4, b = 2, c = 1) ─────────────────────────────────")
let ell = Ellipsoid(4.0, 2.0, 1.0)
    P = hill_tensor(ell, C_iso)
    println()
    print_voigt(P; label="P triaxial")
    # Check symmetries
    err_min = maximum(abs(P[i,j,k,l]-P[j,i,k,l]) for i in 1:3, j in 1:3, k in 1:3, l in 1:3)
    err_maj = maximum(abs(P[i,j,k,l]-P[k,l,i,j]) for i in 1:3, j in 1:3, k in 1:3, l in 1:3)
    @printf "\n  Symmetry: minor err=%.2e,  major err=%.2e\n" err_min err_maj
end

# ═══════════════════════════════════════════════════════════════════════════
println("\n", "="^70)
println("  § 2  ANISOTROPIC MATRIX — residue theorem vs DECUHR")
println("="^70)

# Cubic crystal-like matrix (Zener ratio A = 2C44/(C11−C12) ≠ 1)
# C11=250, C12=100, C44=80 GPa  →  A ≈ 1.067
let
    C11, C12, C44 = 250e3, 100e3, 80e3   # MPa
    C_arr = zeros(3,3,3,3)
    for (I,(i,j)) in enumerate(voigt_idx), (J,(k,l)) in enumerate(voigt_idx)
        v = [C11 C12 C12 0 0 0;
             C12 C11 C12 0 0 0;
             C12 C12 C11 0 0 0;
             0   0   0   C44 0 0;
             0   0   0   0 C44 0;
             0   0   0   0 0 C44][I,J]
        C_arr[i,j,k,l] = C_arr[j,i,k,l] = C_arr[i,j,l,k] = C_arr[j,i,l,k] = v
    end
    C_cubic = Tens(C_arr)

    @printf "\nCubic matrix: C11=%.0f, C12=%.0f, C44=%.0f MPa\n" C11 C12 C44
    @printf "Zener ratio A = 2C44/(C11−C12) = %.4f  (iso = 1.0)\n" 2C44/(C11-C12)

    # Warm-up (avoid counting compilation)
    let _ell = Ellipsoid(2.0, 1.0, 1.0)
        hill_tensor(_ell, C_cubic; method=:residues, abstol=1e-8, reltol=1e-6)
        hill_tensor(_ell, C_cubic; method=:decuhr,  abstol=1e-8, reltol=1e-6)
    end

    ell = Ellipsoid(3.0, 1.0, 1.0)
    println("\n  Prolate spheroid (a=3, b=c=1) in cubic matrix:")

    t_res = @elapsed P_res = hill_tensor(ell, C_cubic; method=:residues, abstol=1e-8, reltol=1e-6)
    t_dcr = @elapsed P_dcr = hill_tensor(ell, C_cubic; method=:decuhr,  abstol=1e-8, reltol=1e-6)

    @printf "  :residues  time=%.4f s\n" t_res
    @printf "  :decuhr   time=%.4f s\n" t_dcr

    max_err = maximum(abs(P_res[i,j,k,l]-P_dcr[i,j,k,l]) for i in 1:3, j in 1:3, k in 1:3, l in 1:3)
    @printf "  Max |P_residue − P_decuhr| = %.3e MPa⁻¹\n" max_err

    println("\n  Voigt[:residues] vs Voigt[:decuhr]  (×10⁻⁶ MPa⁻¹):")
    println("        label     :residues    :decuhr      diff")
    for (I,(i,j)) in enumerate(voigt_idx), (J,(k,l)) in enumerate(voigt_idx)
        J >= I || continue
        vr = P_res[i,j,k,l]; vd = P_dcr[i,j,k,l]
        abs(vr) > 1e-12 || abs(vd) > 1e-12 || continue
        @printf "  P[%s,%s]: %10.4f  %10.4f  %10.2e\n" voigt_lab[I] voigt_lab[J] 1e6*vr 1e6*vd abs(vr-vd)
    end
end

# ═══════════════════════════════════════════════════════════════════════════
println("\n", "="^70)
println("  § 3  2D PLANE STRAIN")
println("="^70)

let
    C_iso2 = TensISO{2}(3k_ref, 2μ_ref)

    println("\n── Circle  (r = 1)  isotropic matrix:")
    ell = Ellipsoid(1.0; dim=2)
    P = hill_tensor(ell, C_iso2)
    @printf "  P[1,1,1,1] = %12.9e,  P[2,2,2,2] = %12.9e  (should be equal)\n" P[1,1,1,1] P[2,2,2,2]
    @printf "  Isotropy error = %.2e\n" abs(P[1,1,1,1]-P[2,2,2,2])

    println("\n── Ellipse  (a=4, b=1)  isotropic matrix:")
    ell = Ellipsoid(4.0, 1.0)
    P = hill_tensor(ell, C_iso2)
    @printf "  P[1,1,1,1] = %12.9e  (major-axis direction)\n" P[1,1,1,1]
    @printf "  P[2,2,2,2] = %12.9e  (minor-axis direction)\n" P[2,2,2,2]
    @printf "  P[1,2,1,2] = %12.9e\n" P[1,2,1,2]

    println("\n── Ellipse  (a=4, b=1)  orthorhombic matrix (E1=100, E2=200, ν12=0.3, G12=40 GPa):")
    E1, E2, ν12, G12 = 100e3, 200e3, 0.3, 40e3
    ν21 = ν12 * E2/E1
    D   = 1 - ν12*ν21
    C2_arr = zeros(2,2,2,2)
    C2_arr[1,1,1,1] = E1/D;  C2_arr[2,2,2,2] = E2/D
    for idx in ((1,1,2,2),(2,2,1,1)); C2_arr[idx...] = ν12*E2/D; end
    for idx in ((1,2,1,2),(1,2,2,1),(2,1,1,2),(2,1,2,1)); C2_arr[idx...] = G12; end
    C2_aniso = Tens(C2_arr)

    ell = Ellipsoid(4.0, 1.0)
    P = hill_tensor(ell, C2_aniso; abstol=1e-8, reltol=1e-6)
    @printf "  P[1,1,1,1] = %12.9e MPa⁻¹\n" P[1,1,1,1]
    @printf "  P[2,2,2,2] = %12.9e MPa⁻¹\n" P[2,2,2,2]
    @printf "  P[1,2,1,2] = %12.9e MPa⁻¹\n" P[1,2,1,2]
end

# ═══════════════════════════════════════════════════════════════════════════
println("\n", "="^70)
println("  § 4  ESHELBY TENSOR  S^E = P : C₀")
println("="^70)

println("\n  Sphere, isotropic matrix — Eshelby (1957) analytical values:")
let ell = Ellipsoid(1.0)
    P = hill_tensor(ell, C_iso)
    S = dcontract(P, C_iso)

    # Analytical formulas
    S1111_th = (7 - 5ν_ref) / (15*(1 - ν_ref))
    S1122_th = (5ν_ref - 1) / (15*(1 - ν_ref))
    S1212_th = (4 - 5ν_ref) / (15*(1 - ν_ref))

    println("\n  Component       Computed       Analytical     Error")
    @printf "  S[1,1,1,1]  %12.8f   %12.8f   %.2e\n" S[1,1,1,1] S1111_th abs(S[1,1,1,1]-S1111_th)
    @printf "  S[1,1,2,2]  %12.8f   %12.8f   %.2e\n" S[1,1,2,2] S1122_th abs(S[1,1,2,2]-S1122_th)
    @printf "  S[1,2,1,2]  %12.8f   %12.8f   %.2e\n" S[1,2,1,2] S1212_th abs(S[1,2,1,2]-S1212_th)

    println("\n  Prolate spheroid  (a=5, b=c=1) — Eshelby tensor:")
    P2 = hill_tensor(Ellipsoid(5.0, 1.0, 1.0), C_iso)
    S2 = dcontract(P2, C_iso)
    @printf "  S[1,1,1,1] = %.7f  (axial)\n" S2[1,1,1,1]
    @printf "  S[2,2,2,2] = %.7f  (transverse)\n" S2[2,2,2,2]
    @printf "  S[1,1,2,2] = %.7f\n" S2[1,1,2,2]
    @printf "  Transverse isotropy: S[2222]=S[3333]? err=%.2e\n" abs(S2[2,2,2,2]-S2[3,3,3,3])
end

# ═══════════════════════════════════════════════════════════════════════════
println("\n", "="^70)
println("  § 5  DILUTE HOMOGENIZATION — spherical pores in an isotropic matrix")
println("="^70)
println("  Dilute estimate:  C_eff = C₀ + f δC : (I + P:δC)⁻¹")
println("  For voids: C_i = 0  →  δC = −C₀,  A_dil = (I − S^E)⁻¹")
println()

let
    f  = 0.05  # 5 % porosity
    ell = Ellipsoid(1.0)
    P  = hill_tensor(ell, C_iso)

    M₀ = mandel(C_iso)                   # 6×6 Mandel of C₀
    Mδ = -M₀                              # δC = −C₀ for voids
    MP = mandel(P)                        # 6×6 Mandel of P

    I6  = Matrix{Float64}(I, 6, 6)
    M_A = inv(I6 + MP * Mδ)              # Mandel localization A = (I + P:δC)⁻¹
    M_eff = M₀ + f * Mδ * M_A            # Dilute effective stiffness (Mandel)
    S_eff = inv(M_eff)                    # Effective compliance (Mandel)

    E1_eff = 1 / S_eff[1,1]
    ν12_eff = -S_eff[1,2] / S_eff[1,1]
    μ_eff  = 1 / S_eff[6,6]

    @printf "  Void fraction f = %.2f\n\n" f
    @printf "  Matrix:    E = %.1f MPa,  ν = %.3f,  μ = %.1f MPa\n" E_ref ν_ref μ_ref
    @printf "  Effective: E = %.1f MPa  (E_eff/E₀ = %.4f)\n" E1_eff E1_eff/E_ref
    @printf "  Effective: ν = %.4f\n" ν12_eff
    @printf "  Effective: μ = %.1f MPa  (μ_eff/μ₀ = %.4f)\n" μ_eff μ_eff/μ_ref

    # Analytical dilute result for spherical pores (Eshelby 1957):
    # k_eff = k₀ (1 − f*(3k₀+4μ₀)/(k₀+4μ₀/3*(1−f))) ≈ k₀(1−f*β_k)
    # μ_eff = μ₀ (1 − 5f*(3k₀+4μ₀)/(9k₀+8μ₀)) + O(f²)  (dilute limit)
    k₀ = k_ref;  μ₀ = μ_ref
    μ_eff_th = μ₀ * (1 - 5f*(3k₀+4μ₀)/(9k₀+8μ₀))
    @printf "\n  Analytical μ_eff (dilute) = %.1f MPa  (err=%.2e)\n" μ_eff_th abs(μ_eff-μ_eff_th)
end

println()
println("="^70)
