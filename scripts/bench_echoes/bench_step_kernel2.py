"""
bench_step_kernel2.py — extract Python ECHOES' per-layer trapezoidal
Volterra matrices (k, mu) for the step-activated layered sphere setup.
We use the `visco_mat(name, layer)` accessor to read the internal
storage that the recurrence operates on.
"""
from numpy import *
from echoes import *
import json
set_printoptions(precision=8, suppress=True, linewidth=200)

E0 = 1.0; nu0 = 0.2; C0 = stiff_Enu(E0, nu0)
f0 = 0.6; eta0 = 0.2; gamma0 = 0.133
E1 = 5.0; nu1 = 0.3; C1 = stiff_Enu(E1, nu1)
finf = 0.3; eta1 = 1.0; gamma1 = 1.67
Ep = E0 * 1.0e-8; nup = 0.2; Cp = stiff_Enu(Ep, nup); fp = 1 - f0 - finf

R0 = lambda t, tp: 3*C0.k*exp(-(t-tp)/eta0)*J4 + 2*C0.mu*exp(-(t-tp)/gamma0)*K4
R1 = lambda t, tp: 3*C1.k*exp(-(t-tp)/eta1)*J4 + 2*C1.mu*exp(-(t-tp)/gamma1)*K4

# Solidification with N=3, alpha=4.
N = 3
alpha = 4.0
F = array([(i+0.5)*finf/N for i in range(N)])
lT = array([(f/(finf-f))**(1.0/alpha) for f in F])

T = array([0.5, 0.7, 1.0, 1.5, 2.0])
nT = len(T)
t0 = T[0]

ver = rve(matrix="MATRIX")
ver["MATRIX"] = ellipsoid(spherical, fraction=f0,
                          prop={"C": C0},
                          visco_prop={"C": (R0, RELAXATION)})

ver["INCLUSION"] = sphere_nlayers(
    radius=1.0,
    layer_fractions=[fp] + [finf/N for _ in range(N)],
    fraction=finf+fp,
    prop={"C": [Cp] + [C1 if t0>=lT[N-1-i] else Cp for i in range(N)]},
    visco_prop={"C": [(lambda _,__: Cp.array, RELAXATION)] +
        [(lambda t, tp, lt=lT[N-1-i]: R1(t, tp) if tp>=lt else Cp.array, RELAXATION)
         for i in range(N)]},
)

# Trigger build.
homogenize_visco(prop="C", rve=ver, time_series=T, scheme=MT,
                 maxnb=100, epsrel=1.e-8, verbose=False)

# Per-layer (k, mu) trapezoidal matrices, Python-indexed 0..N.
out = {"T": T.tolist(), "lT": lT.tolist(), "t0": t0, "N": N,
       "layers_k": [], "layers_mu": []}

# Layer 0 = innermost (pore), layer N = outermost (lT[0] = earliest activation).
for k in range(N+1):
    Mk = ver["INCLUSION"].visco_mat("Ck", k)
    Mmu = ver["INCLUSION"].visco_mat("Cmu", k)
    print(f"\n--- Python ECHOES layer {k} ---")
    print(f"k Volterra (n×n):"); print(Mk)
    print(f"mu Volterra (n×n):"); print(Mmu)
    out["layers_k"].append(Mk.tolist())
    out["layers_mu"].append(Mmu.tolist())

import os
out_path = os.path.join(os.path.dirname(__file__), "bench_step_layer_kernel_python.json")
with open(out_path, "w") as f:
    json.dump(out, f, indent=2)
print(f"\nSaved to {out_path}")
