module MeanFieldHomogenizationSymPyExt

using MeanFieldHomogenization
using SymPy

# ──────────────────────────────────────────────────────────────────────────────
#  Symbolic closed forms via SymPy's `elliptic_{k,e,f}`
#
#  SymPy conventions (identical to `MeanFieldHomogenization.Elliptic` ≡ `Elliptic.jl`):
#  the parameter is `m = k²` (not the modulus `k`).
#
#  Without this extension, the generic AGM path would unfold ~60 nested
#  `sqrt` expressions on a `Sym` input and overwhelm SymPy's pretty-printer.
# ──────────────────────────────────────────────────────────────────────────────

MeanFieldHomogenization.Elliptic.ell_K(m::Sym) = sympy.elliptic_k(m)
MeanFieldHomogenization.Elliptic.ell_E(m::Sym) = sympy.elliptic_e(m)

MeanFieldHomogenization.Elliptic.ell_F(φ::Sym, m::Sym) = sympy.elliptic_f(φ, m)
MeanFieldHomogenization.Elliptic.ell_E(φ::Sym, m::Sym) = sympy.elliptic_e(φ, m)

# Mixed-type cases — promote the non-Sym argument
MeanFieldHomogenization.Elliptic.ell_F(φ::Sym, m::Number) = sympy.elliptic_f(φ, Sym(m))
MeanFieldHomogenization.Elliptic.ell_F(φ::Number, m::Sym) = sympy.elliptic_f(Sym(φ), m)
MeanFieldHomogenization.Elliptic.ell_E(φ::Sym, m::Number) = sympy.elliptic_e(φ, Sym(m))
MeanFieldHomogenization.Elliptic.ell_E(φ::Number, m::Sym) = sympy.elliptic_e(Sym(φ), m)

end
