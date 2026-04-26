"""
bench_layered_alv_step.py — cross-check layered-sphere ALV with the
step-activation pattern of `fluage_echoes_solid.py`.

Setup mirrors script 37 but with a small N (3 layers) and 5 time steps
so the (n×n) Volterra matrices are printable.  At t₀ = 0.5, layer 1
(innermost) is already activated (t_set ≈ 0.46), layer 2 not yet
(t_set ≈ 0.79).
"""

from numpy import *
from echoes import *
import json
set_printoptions(precision=8, suppress=True, linewidth=200)

# ── Material setup matching script 37 ─────────────────────────────────────

E0 = 1.0; nu0 = 0.2; C0 = stiff_Enu(E0, nu0)
f0 = 0.6; eta0 = 0.2; gamma0 = 0.133
E1 = 5.0; nu1 = 0.3; C1 = stiff_Enu(E1, nu1)
finf = 0.3; eta1 = 1.0; gamma1 = 1.67
Ep = E0 * 1.0e-8; nup = 0.2; Cp = stiff_Enu(Ep, nup); fp = 1 - f0 - finf

R0 = lambda t, tp: 3*C0.k*exp(-(t-tp)/eta0)*J4 + 2*C0.mu*exp(-(t-tp)/gamma0)*K4
R1 = lambda t, tp: 3*C1.k*exp(-(t-tp)/eta1)*J4 + 2*C1.mu*exp(-(t-tp)/gamma1)*K4

# Solidification kinetics
N = 3
alpha = 4.0
F = array([(i+0.5)*finf/N for i in range(N)])
lT = array([(f/(finf-f))**(1.0/alpha) for f in F])
print("setting times lT =", lT)

# Time grid (small).
T = array([0.5, 0.7, 1.0, 1.5, 2.0])
nT = len(T)
t0 = T[0]
print("t0 =", t0)

# Build layered RVE matching script 37 :layers.
ver = rve(matrix="MATRIX")
ver["MATRIX"] = ellipsoid(spherical, fraction=f0,
                          prop={"C": C0},
                          visco_prop={"C": (R0, RELAXATION)})

# Layers: pore + N solidifying shells.
# Outer-most is layer N (innermost solidifying = lT[0]); innermost is
# pore.  Lambda capture trick: lt=lT[N-1-i].
ver["INCLUSION"] = sphere_nlayers(
    radius=1.0,
    layer_fractions=[fp] + [finf/N for _ in range(N)],
    fraction=finf+fp,
    prop={"C": [Cp] + [C1 if t0>=lT[N-1-i] else Cp for i in range(N)]},
    visco_prop={"C": [(lambda _,__: Cp.array, RELAXATION)] +
        [(lambda t, tp, lt=lT[N-1-i]: R1(t, tp) if tp>=lt else Cp.array, RELAXATION)
         for i in range(N)]},
)

# Trigger MT.
V = homogenize_visco(prop="C", rve=ver, time_series=T, scheme=MT,
                     maxnb=100, epsrel=1.e-8, verbose=False)

print("\n=== Per-layer alpha (3K) and beta (2mu) Volterra n×n matrices ===")
out = {"T": T.tolist(), "lT": lT.tolist(), "t0": t0,
       "k0": C0.k, "mu0": C0.mu, "eta0": eta0, "gamma0": gamma0,
       "k1": C1.k, "mu1": C1.mu, "eta1": eta1, "gamma1": gamma1,
       "kp": Cp.k, "mup": Cp.mu, "fp": fp, "finf": finf, "f0": f0,
       "N": N, "alpha": alpha,
       "layers_alpha": [], "layers_beta": []}
for k in range(N+1):
    raw = ver["INCLUSION"].layer_visco_eE(k)
    a, b = visco_paramsym(raw, ISO)
    print(f"\n--- Layer {k} ---")
    print("alpha:")
    print(a)
    print("beta:")
    print(b)
    out["layers_alpha"].append(a.tolist())
    out["layers_beta"].append(b.tolist())

# Effective creep compliance.
J = linalg.inv(V)
JE = sum(J[::6, ::6], 1)
print("\nEffective creep compliance J^E_eff(t, t0):")
print(JE)
out["JE_eff"] = JE.tolist()

import os
out_path = os.path.join(os.path.dirname(__file__), "bench_layered_alv_step_python.json")
with open(out_path, "w") as f:
    json.dump(out, f, indent=2)
print(f"\nSaved Python step reference to {out_path}")
