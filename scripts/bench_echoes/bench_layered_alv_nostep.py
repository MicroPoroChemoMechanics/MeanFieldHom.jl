"""
bench_layered_alv_nostep.py — same geometry as the step benchmark
(N=3 solidifying + 1 pore = 4 layers) but ALL layers elastic (no
step), so any divergence must come from layered-sphere recurrence
itself with multi-layer (>2) configurations.
"""
from numpy import *
from echoes import *
import json
set_printoptions(precision=8, suppress=True, linewidth=200)

E0 = 1.0; nu0 = 0.2; C0 = stiff_Enu(E0, nu0)
f0 = 0.6; eta0 = 0.2; gamma0 = 0.133
E1 = 5.0; nu1 = 0.3; C1 = stiff_Enu(E1, nu1)
finf = 0.3
Ep = E0 * 1.0e-8; nup = 0.2; Cp = stiff_Enu(Ep, nup); fp = 1 - f0 - finf

R0 = lambda t, tp: 3*C0.k*exp(-(t-tp)/eta0)*J4 + 2*C0.mu*exp(-(t-tp)/gamma0)*K4

N = 3
T = array([0.5, 0.7, 1.0, 1.5, 2.0])
nT = len(T)

ver = rve(matrix="MATRIX")
ver["MATRIX"] = ellipsoid(spherical, fraction=f0,
                          prop={"C": C0},
                          visco_prop={"C": (R0, RELAXATION)})

# ALL layers elastic, no step.  Innermost = pore, outermost three = R1 stiff.
ver["INCLUSION"] = sphere_nlayers(
    radius=1.0,
    layer_fractions=[fp] + [finf/N for _ in range(N)],
    fraction=finf+fp,
    prop={"C": [Cp, C1, C1, C1]},
)
homogenize_visco(prop="C", rve=ver, time_series=T, scheme=MT,
                 maxnb=100, epsrel=1.e-8, verbose=False)

out = {"T": T.tolist(),
       "k0": C0.k, "mu0": C0.mu, "eta0": eta0, "gamma0": gamma0,
       "k1": C1.k, "mu1": C1.mu,
       "kp": Cp.k, "mup": Cp.mu, "fp": fp, "finf": finf, "f0": f0,
       "N": N,
       "layers_alpha": [], "layers_beta": []}
for k in range(N+1):
    raw = ver["INCLUSION"].layer_visco_eE(k)
    a, b = visco_paramsym(raw, ISO)
    print(f"\n--- Layer {k} ---")
    print("alpha (3K):")
    print(a)
    print("beta (2mu):")
    print(b)
    out["layers_alpha"].append(a.tolist())
    out["layers_beta"].append(b.tolist())

import os
out_path = os.path.join(os.path.dirname(__file__), "bench_layered_alv_nostep_python.json")
with open(out_path, "w") as f:
    json.dump(out, f, indent=2)
print("Saved", out_path)
