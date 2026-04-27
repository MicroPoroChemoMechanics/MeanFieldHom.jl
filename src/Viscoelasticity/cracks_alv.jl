# =============================================================================
#  cracks_alv.jl — ageing linear viscoelastic cracks (skeleton + roadmap).
#
#  ECHOES C++ exposes ALV cracks via the `crack(shape, interf_visco_prop, …)`
#  Python interface, see
#  `c:/Users/jf.barthelemy/VSCode_workspace/Echoes/echoes_cpp/tests/python/creep/fluage_echoes_cracks.py`.
#
#  The supported geometries / matrices in ECHOES are:
#    * iso ALV matrix `R(t,t')` (relaxation 4-tensor)
#    * thin spheroidal crack with two **interface stiffness** ALV laws
#      `(Rn(t,t'), Rt(t,t'))` (normal and tangential interface compliances)
#
#  The Julia equivalents to develop here are:
#
#    cod_kernel_alv(crack, C₀_law, times)               -> (6n × 6n) Matrix
#    compliance_contribution_alv(crack, C₀_law, times)  -> (6n × 6n) Matrix
#
#  with both routed through the existing scheme dispatcher in
#  `homogenize_alv` via the `CrackDensity` amount on the inclusion phase.
#
#  ── Time–space decoupling for a pure penny crack in an iso ALV matrix ─────
#
#  The elastic COD tensor of a penny crack `(η = 1)` in iso `(K, μ)` is
#  diagonal in the crack basis (l̂, m̂, n̂):
#       B_nn = 16 (1−ν²) / (3π E)
#       B_ll = B_mm = 32 (1−ν²) / (3π E (2−ν))
#
#  In the Volterra n×n algebra with α = 3K, β = 2μ:
#       1 − ν²  =  9K (3K + 4μ) / (4 (3K + μ)²) = α (α + 2β) / (4 (α + β/2)²)
#       E       =  9 K μ / (3 K + μ) = α β / (2 (α + β/2))
#       (1−ν²)/E = (3K + 4μ) / (4μ (3K + μ)) = (α + 2β) / (2β (α + β/2))
#
#  giving
#       B̃_nn = (8 / (3π))  · (α + 2β) ∘ (β ∘ (α + β/2))^{-vol}
#       B̃_ll = (16 / (3π)) · (α + 2β) ∘ (β ∘ (α + β/2) · (2 - ν))^{-vol}   (with 2-ν tensorial)
#
#  Implementing this requires:
#    1. extracting (α, β) from the iso ALV matrix block matrix
#    2. evaluating ν as a Volterra rational fraction in (α, β) → an extra
#       n×n matrix for the (2 - ν) factor in B_ll / B_mm
#    3. assembling the 6×6 block in the crack basis (l̂, m̂, n̂) and rotating
#       to the global Mandel frame at every (t, t').
#
#  Step (3) reduces the problem to four scalar Volterra n×n matrices —
#  one for the diagonal in the crack basis (B_ll, B_mm, B_nn distinct) and
#  one rotation tensor.  The result is then a TI block matrix in the
#  crack normal axis, fitting the existing TI ALV fast path.
#
#  ── Interface-stiffness ALV cracks ──────────────────────────────────────
#
#  When the crack carries `(Rn(t,t'), Rt(t,t'))` interface stiffnesses
#  (cf. `fluage_echoes_cracks.py` setup), the COD tensor at every (t,t')
#  is:
#       B = (4η / 3) · diag(1/Rt, 1/Rt, 1/Rn)   (in the crack basis)
#
#  whose Volterra-discretised version is straightforward:
#       B̃_n = (4η / 3) · Rn^{-vol}   (n × n scalar Volterra inverse)
#       B̃_t = (4η / 3) · Rt^{-vol}
#  Two scalar Volterra inverses, then assembly into a TI 6n×6n block
#  matrix in the crack normal axis.
#
#  ── Compliance contribution H̃ ──────────────────────────────────────────
#
#  Following the elastic bridge:
#       H̃ = (3/4) · n̂ ⊗ˢ B̃ ⊗ˢ n̂   (elliptic / penny)
#       H̃ = (2/π) · n̂ ⊗ˢ B̃ ⊗ˢ n̂   (ribbon)
#
#  Symbolic-tensor product on each (i, j) Volterra entry separately —
#  identical to the elastic case at each time pair.
#
#  ── Integration into homogenize_alv ─────────────────────────────────────
#
#  Add `_inclusion_alv_quantities(crack::AbstractCrack, ...)` that returns
#  `(C_r ≡ 0, A_dut ≡ 𝟙, N_dut ≡ H̃, P_r)` for cracks (the contribution to
#  effective compliance), reusing the existing CrackDensity dispatch from
#  the elastic Schemes module.
#
#  The Mori-Tanaka / SC / PCW algebra in `schemes_alv.jl` is unchanged —
#  cracks just plug in as one more inhomogeneity with H̃ as their
#  compliance contribution.
#
#  ── Status ──────────────────────────────────────────────────────────────
#
#  * Pure penny crack in iso ALV matrix: implementation pending.
#  * Interface-stiffness crack (Rn, Rt): implementation pending.
#  * Test against `fluage_echoes_cracks.py` benchmark: pending.
# =============================================================================

# Placeholder exports — concrete definitions to be added in the follow-up.
# `cod_kernel_alv` / `compliance_contribution_alv` will mirror
# `cod_tensor` / `compliance_contribution` from `Cracks/`.
