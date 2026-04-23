module MeanFieldHomSymPyExt

using MeanFieldHom
using SymPy

# ──────────────────────────────────────────────────────────────────────────────
#  Symbolic closed forms via SymPy's `elliptic_{k,e,f}`
#
#  SymPy conventions (identical to `MeanFieldHom.Elliptic` ≡ `Elliptic.jl`):
#  the parameter is `m = k²` (not the modulus `k`).
#
#  Without this extension, the generic AGM path would unfold ~60 nested
#  `sqrt` expressions on a `Sym` input and overwhelm SymPy's pretty-printer.
# ──────────────────────────────────────────────────────────────────────────────

MeanFieldHom.Elliptic.ell_K(m::Sym) = sympy.elliptic_k(m)
MeanFieldHom.Elliptic.ell_E(m::Sym) = sympy.elliptic_e(m)

MeanFieldHom.Elliptic.ell_F(φ::Sym, m::Sym) = sympy.elliptic_f(φ, m)
MeanFieldHom.Elliptic.ell_E(φ::Sym, m::Sym) = sympy.elliptic_e(φ, m)

# Mixed-type cases — promote the non-Sym argument
MeanFieldHom.Elliptic.ell_F(φ::Sym, m::Number) = sympy.elliptic_f(φ, Sym(m))
MeanFieldHom.Elliptic.ell_F(φ::Number, m::Sym) = sympy.elliptic_f(Sym(φ), m)
MeanFieldHom.Elliptic.ell_E(φ::Sym, m::Number) = sympy.elliptic_e(φ, Sym(m))
MeanFieldHom.Elliptic.ell_E(φ::Number, m::Sym) = sympy.elliptic_e(Sym(φ), m)

end
