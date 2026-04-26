"""N=2 step layers + Maxwell matrix, no pore."""
from numpy import *
from echoes import *
import json
set_printoptions(precision=8, suppress=True, linewidth=200)

E0 = 1.0; nu0 = 0.2; C0 = stiff_Enu(E0, nu0)
eta0 = 0.2; gamma0 = 0.133
E1 = 5.0; nu1 = 0.3; C1 = stiff_Enu(E1, nu1)
eta1 = 1.0; gamma1 = 1.67
Ep = 1e-8; nup = 0.2; Cp = stiff_Enu(Ep, nup)

R0 = lambda t, tp: 3*C0.k*exp(-(t-tp)/eta0)*J4 + 2*C0.mu*exp(-(t-tp)/gamma0)*K4
R1 = lambda t, tp: 3*C1.k*exp(-(t-tp)/eta1)*J4 + 2*C1.mu*exp(-(t-tp)/gamma1)*K4

T = array([0.5, 0.7, 1.0, 1.5, 2.0])
t_set_inner = 1.0    # layer 1 (inner) activates at t=1.0
t_set_outer = 0.7    # layer 2 (outer) activates at t=0.7

ver = rve(matrix="MAT")
ver["MAT"] = ellipsoid(spherical, fraction=0.6,
                       prop={"C": C0}, visco_prop={"C": (R0, RELAXATION)})
ver["INC"] = sphere_nlayers(radius=1.0,
    layer_fractions=[0.5, 0.5],
    fraction=0.4,
    prop={"C": [Cp, Cp]},   # both pore at t=0.5 (first time)
    visco_prop={"C": [
        (lambda t, tp: R1(t, tp) if tp >= t_set_inner else Cp.array, RELAXATION),
        (lambda t, tp: R1(t, tp) if tp >= t_set_outer else Cp.array, RELAXATION),
    ]})
homogenize_visco(prop="C", rve=ver, time_series=T, scheme=MT, maxnb=100, epsrel=1e-8, verbose=False)

out = {"T": T.tolist(), "k0": C0.k, "mu0": C0.mu, "eta0": eta0, "gamma0": gamma0,
       "k1": C1.k, "mu1": C1.mu, "eta1": eta1, "gamma1": gamma1,
       "kp": Cp.k, "mup": Cp.mu, "t_set_inner": t_set_inner, "t_set_outer": t_set_outer,
       "layers_alpha": [], "layers_beta": []}
for k in range(2):
    raw = ver["INC"].layer_visco_eE(k)
    a, b = visco_paramsym(raw, ISO)
    print(f"\n--- Layer {k} ---")
    print("alpha (3K):"); print(a)
    print("beta (2mu):"); print(b)
    out["layers_alpha"].append(a.tolist())
    out["layers_beta"].append(b.tolist())

import os
out_path = os.path.join(os.path.dirname(__file__), "bench_step_n2_python.json")
with open(out_path, "w") as f:
    json.dump(out, f, indent=2)
print("Saved", out_path)
