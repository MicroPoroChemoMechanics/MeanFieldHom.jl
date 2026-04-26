"""Sweep N=2,3,4 elastic + Maxwell matrix to find where it breaks."""
from numpy import *
from echoes import *
import json
set_printoptions(precision=8, suppress=True, linewidth=200)

E0 = 1.0; nu0 = 0.2; C0 = stiff_Enu(E0, nu0)
eta0 = 0.2; gamma0 = 0.133
R0 = lambda t, tp: 3*C0.k*exp(-(t-tp)/eta0)*J4 + 2*C0.mu*exp(-(t-tp)/gamma0)*K4
T = array([0.5, 0.7, 1.0, 1.5, 2.0])

ks_all  = [0.5, 1.0, 2.0, 3.0]
mus_all = [0.3, 0.6, 1.0, 1.5]

results = {"T": T.tolist(), "k0": C0.k, "mu0": C0.mu, "eta0": eta0, "gamma0": gamma0, "Ns": []}

for N in [2, 3, 4]:
    Cs = [stiff_kmu(ks_all[i], mus_all[i]) for i in range(N)]
    finf = 0.4
    ver = rve(matrix="MAT")
    ver["MAT"] = ellipsoid(spherical, fraction=1.0 - finf,
                           prop={"C": C0}, visco_prop={"C": (R0, RELAXATION)})
    ver["INC"] = sphere_nlayers(radius=1.0,
        layer_fractions=[1.0/N for _ in range(N)],
        fraction=finf, prop={"C": Cs})
    homogenize_visco(prop="C", rve=ver, time_series=T, scheme=MT, maxnb=100, epsrel=1e-8, verbose=False)
    item = {"N": N, "ks": ks_all[:N], "mus": mus_all[:N], "alpha": [], "beta": []}
    for k in range(N):
        raw = ver["INC"].layer_visco_eE(k)
        a, b = visco_paramsym(raw, ISO)
        item["alpha"].append(a.tolist())
        item["beta"].append(b.tolist())
    results["Ns"].append(item)

import os
out_path = os.path.join(os.path.dirname(__file__), "debug_n_sweep_python.json")
with open(out_path, "w") as f:
    json.dump(results, f, indent=2)
print("Saved", out_path)
