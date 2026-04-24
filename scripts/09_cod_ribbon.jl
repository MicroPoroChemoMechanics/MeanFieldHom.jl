# =============================================================================
#  scripts/04_cod_ribbon.jl
#
#  Compare ribbon-crack COD tensors obtained from :analytical (ISO),
#  :residues (ISO promoted to Tens) and :decuhr (idem), and verify the
#  consistency relation B^R = (3π/8) lim_{η→0} B^ℰ(η) numerically.
# =============================================================================

import Pkg
Pkg.activate(joinpath(@__DIR__, ".."); io=devnull)

using MeanFieldHom
using TensND
using LinearAlgebra

E, ν = 1.0, 0.25
k, μ = E / (3(1-2ν)), E / (2(1+ν))
C_iso = TensISO{3}(3k, 2μ)

# Promote to generic Tens to force the numerical dispatch
C_generic = zeros(3,3,3,3)
for i=1:3, j=1:3, p=1:3, q=1:3
    C_generic[i,j,p,q] = C_iso[i,j,p,q]
end
C_gen = Tens(C_generic, CanonicalBasis{3,Float64}())

r = RibbonCrack(1.0)
println("Analytical (TensISO): ", cod_tensor(r, C_iso))
println("Residue    (generic): ", cod_tensor(r, C_gen; method=:residues))
println("DECUHR     (generic): ", cod_tensor(r, C_gen; method=:decuhr))

# Consistency check: (3π/8)·B^E(η→0) ≈ B^R   on nn and mm components
B_r = cod_tensor(r, C_iso)
println("\nRibbon reference  B^R[3,3] = $(B_r[3,3])")
for η in [1e-2, 1e-3, 1e-4]
    ec = EllipticCrack(1.0, η)
    B_e = cod_tensor(ec, C_iso)
    println("(3π/8) · B^E(η=$η)[3,3] = $((3π/8) * B_e[3,3])")
end
