"""
bench_layered_alv.py — minimal cross-check of per-layer ALV concentration
tensors between the reference ECHOES Python implementation and the Julia
`MeanFieldHom.jl` re-implementation.

We pin a small, fully-explicit setup (Maxwell matrix + 2 elastic layers
+ 4 time points) so the resulting (6n x 6n) Volterra matrices fit in
the terminal. The output JSON file is consumed by
`bench_layered_alv.jl` for side-by-side comparison.
"""

from numpy import *
from echoes import *
import json
set_printoptions(precision=8, suppress=True)

# ── Material setup ─────────────────────────────────────────────────────────
# Matrix: Maxwell relaxation R(t,t') = 3 k0 exp(-(t-t')/eta_k) J + 2 mu0 exp(-(t-t')/eta_mu) K
k0 = 1.0; mu0 = 0.5
eta_k0 = 0.6; eta_mu0 = 2.0
R0_lam = lambda t, tp: (3.0 * k0 * exp(-(t - tp) / eta_k0)) * J4 \
                        + (2.0 * mu0 * exp(-(t - tp) / eta_mu0)) * K4

# Two elastic inclusion layers (Heaviside relaxation kernels).
k1 = 2.0; mu1 = 1.0      # core
k2 = 3.0; mu2 = 1.5      # shell
C1 = stiff_kmu(k1, mu1)
C2 = stiff_kmu(k2, mu2)
R1 = lambda t, tp: (3.0 * k1) * J4 + (2.0 * mu1) * K4 if t >= tp else zeros((3, 3, 3, 3))
R2 = lambda t, tp: (3.0 * k2) * J4 + (2.0 * mu2) * K4 if t >= tp else zeros((3, 3, 3, 3))

# Geometry: composite sphere with two layers of equal volume fraction
# (radii at (0.5)^(1/3) and 1.0).
fp = 0.0
finf = 0.4               # total inclusion volume fraction in the matrix
N_layers = 2
layer_fracs = [0.5, 0.5]    # within the sphere

# Time grid (small so we can print everything).
T = array([0.0, 0.5, 1.0, 1.5, 2.0])
nT = len(T)

# Build RVE.
ver = rve(matrix="MATRIX")
ver["MATRIX"] = ellipsoid(spherical, fraction=1.0 - finf,
                          prop={"C": tensor(R0_lam(0.0, 0.0))},
                          visco_prop={"C": (R0_lam, RELAXATION)})

ver["INCLUSION"] = sphere_nlayers(
    radius=1.0,
    layer_fractions=layer_fracs,
    fraction=finf,
    prop={"C": [C1, C2]},
    visco_prop={"C": [(R1, RELAXATION), (R2, RELAXATION)]},
)

# Trigger viscoelastic homogenization (this populates the internal
# Volterra matrices).
V = homogenize_visco(prop="C", rve=ver, time_series=T, scheme=MT,
                     maxnb=100, epsrel=1.e-8, verbose=False)

# Extract per-layer raw (dilute-level) concentration tensors.
# Each is a (6n x 6n) ndarray.
layer0 = ver["INCLUSION"].layer_visco_eE(0)
layer1 = ver["INCLUSION"].layer_visco_eE(1)
matrix_eE = ver["MATRIX"].visco_eE

print("layer_visco_eE shape:", layer0.shape)
print("Time grid:", T)

# Iso (alpha, beta) extraction via visco_paramsym(_, ISO):
# visco_paramsym returns (3K_volterra, 2mu_volterra) per the ECHOES convention.
def iso_alpha_beta(M):
    a3K, b2mu = visco_paramsym(M, ISO)
    # Convert to our (alpha = M[1,1] + 2 M[1,2] / no, wait — this is the
    # 'iso parameter' API which already returns scalar Volterra blocks).
    return a3K, b2mu

alpha0, beta0 = iso_alpha_beta(layer0)
alpha1, beta1 = iso_alpha_beta(layer1)

print("\n=== Layer 0 (core, k=2, mu=1) ===")
print("alpha (3K volterra):"); print(alpha0)
print("beta (2mu volterra):"); print(beta0)
print("\n=== Layer 1 (shell, k=3, mu=1.5) ===")
print("alpha (3K volterra):"); print(alpha1)
print("beta (2mu volterra):"); print(beta1)

# Save to JSON for Julia consumption.
out = {
    "T": T.tolist(),
    "k0": k0, "mu0": mu0, "eta_k0": eta_k0, "eta_mu0": eta_mu0,
    "k1": k1, "mu1": mu1, "k2": k2, "mu2": mu2,
    "fp": fp, "finf": finf, "layer_fracs": layer_fracs,
    "layer0_eE": layer0.tolist(),
    "layer1_eE": layer1.tolist(),
    "matrix_eE": matrix_eE.tolist(),
    "alpha0": alpha0.tolist(), "beta0": beta0.tolist(),
    "alpha1": alpha1.tolist(), "beta1": beta1.tolist(),
}

import os
out_path = os.path.join(os.path.dirname(__file__), "bench_layered_alv_python.json")
with open(out_path, "w") as f:
    json.dump(out, f, indent=2)
print(f"\nSaved Python reference to {out_path}")
