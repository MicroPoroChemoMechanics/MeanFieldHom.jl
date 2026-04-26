"""
bench_step_kernel.py — dump the trapezoidal Volterra matrix Python
ECHOES builds for a *step-activated* layer kernel.  We want to see
whether ECHOES' build_mat() agrees with Julia's `trapezoidal_matrix`
for an iso step law that switches from Cp to R1 at t_set.
"""
from numpy import *
from echoes import *
import json
set_printoptions(precision=8, suppress=True, linewidth=200)

E0 = 1.0; nu0 = 0.2; C0 = stiff_Enu(E0, nu0)
E1 = 5.0; nu1 = 0.3; C1 = stiff_Enu(E1, nu1)
finf = 0.3; eta1 = 1.0; gamma1 = 1.67
Ep = E0 * 1.0e-8; nup = 0.2; Cp = stiff_Enu(Ep, nup)
fp = 1 - 0.6 - finf

R1_kernel = lambda t, tp: 3*C1.k*exp(-(t-tp)/eta1)*J4 + 2*C1.mu*exp(-(t-tp)/gamma1)*K4

T = array([0.5, 0.7, 1.0, 1.5, 2.0])
nT = len(T)
t_set = 0.669    # falls between T[0]=0.5 and T[1]=0.7

# Build a step law kernel: pore for tp < t_set, R1 for tp >= t_set.
step_kernel = lambda t, tp: R1_kernel(t, tp) if tp >= t_set else Cp.array

# Inject as a "matrix" property in a trivial RVE so ECHOES builds the
# trapezoidal Volterra matrix and we can extract it.
ver = rve(matrix="MAT")
ver["MAT"] = ellipsoid(spherical, fraction=1.0,
                       prop={"C": Cp},   # placeholder
                       visco_prop={"C": (step_kernel, RELAXATION)})

# Trigger build via homogenize_visco (degenerate, but it builds visco_eE
# for the matrix which is the trapezoidal matrix).
homogenize_visco(prop="C", rve=ver, time_series=T, scheme=MT,
                 maxnb=10, epsrel=1e-6, verbose=False)

mat_R = ver["MAT"].visco_eE   # 6n × 6n (relaxation modulus)
print("ECHOES trapezoidal R for step law (6n × 6n):")
print(mat_R)

# Iso decomposition: visco_paramsym(_, ISO) → (3K, 2mu) Volterra n×n.
a3K, b2mu = visco_paramsym(mat_R, ISO)
print("\n=== alpha (3K) Volterra n×n ===")
print(a3K)
print("\n=== beta (2mu) Volterra n×n ===")
print(b2mu)

import os
out_path = os.path.join(os.path.dirname(__file__), "bench_step_kernel_python.json")
with open(out_path, "w") as f:
    json.dump({"T": T.tolist(), "t_set": t_set,
               "alpha": a3K.tolist(), "beta": b2mu.tolist(),
               "k1": C1.k, "mu1": C1.mu, "eta1": eta1, "gamma1": gamma1,
               "kp": Cp.k, "mup": Cp.mu}, f, indent=2)
print(f"\nSaved to {out_path}")
