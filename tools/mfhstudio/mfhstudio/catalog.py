"""What the interface offers, and where each part comes from.

The catalog has two halves, and keeping them apart is what lets the interface
stay usable when Julia is not.

*Form definitions* — which fields a spheroid needs, what a layer looks like,
which lenses exist — are user-interface concerns. They live here, in Python,
and are available immediately.

*Introspected facts* — the list of schemes and the solver options each one
actually reads — are properties of the installed MeanFieldHom and can only come
from the sidecar. They are merged in when it is ready.

The fallback scheme list below exists so the forms render before Julia has
finished loading; it is replaced wholesale by the introspected one, never
merged with it, so a scheme that MFH drops disappears from the interface.
"""

from __future__ import annotations

GEOMETRIES = [
    {
        "name": "Spheroid", "kind": "spheroid", "dim": 3,
        "doc": "Axisymmetric ellipsoid; ω < 1 oblate, ω > 1 prolate, ω = 1 sphere.",
        "fields": [{"name": "omega", "label": "ω (aspect ratio)", "type": "number", "default": 1.0}],
        "angles": 2,
    },
    {
        "name": "Ellipsoid", "kind": "ellipsoid", "dim": 3,
        "doc": "General ellipsoid of semi-axes (a, b, c).",
        "fields": [
            {"name": "a", "label": "a", "type": "number", "default": 1.0},
            {"name": "b", "label": "b", "type": "number", "default": 1.0},
            {"name": "c", "label": "c", "type": "number", "default": 1.0},
        ],
        "angles": 3,
    },
    {
        "name": "Cylinder", "kind": "cylinder", "dim": 3,
        "doc": "Infinite elliptic cylinder of transverse semi-axes (b, c).",
        "fields": [
            {"name": "b", "label": "b", "type": "number", "default": 1.0},
            {"name": "c", "label": "c", "type": "number", "default": 1.0},
        ],
        "angles": 3,
    },
    {
        "name": "PennyCrack", "kind": "penny_crack", "dim": 3,
        "doc": "Circular flat crack of radius a. Enters the RVE with a crack density.",
        "fields": [{"name": "a", "label": "a (radius)", "type": "number", "default": 1.0}],
        "angles": 2, "amount": "density",
    },
    {
        "name": "EllipticCrack", "kind": "elliptic_crack", "dim": 3,
        "doc": "Elliptic flat crack, semi-axes a ≥ b.",
        "fields": [
            {"name": "a", "label": "a", "type": "number", "default": 1.0},
            {"name": "b", "label": "b", "type": "number", "default": 0.5},
        ],
        "angles": 3, "amount": "density",
    },
    {
        "name": "RibbonCrack", "kind": "ribbon_crack", "dim": 3,
        "doc": "Tunnel crack of half-width b, unbounded along its length.",
        "fields": [{"name": "b", "label": "b (half-width)", "type": "number", "default": 1.0}],
        "angles": 3, "amount": "density",
    },
    {
        "name": "LayeredSphere", "kind": "layered_sphere", "dim": 3,
        "doc": "Concentric layers, ascending radii, r = 0 implicit at the center.",
        "fields": [], "layered": True, "angles": 0,
    },
    {
        "name": "LayeredSpheroid", "kind": "layered_spheroid", "dim": 3,
        "doc": "Confocal spheroidal layers. Conduction only.",
        "fields": [
            {"name": "omega", "label": "ω (aspect ratio)", "type": "number", "default": 0.5},
            {"name": "Nseries", "label": "N (series order)", "type": "integer", "default": 5},
        ],
        "layered": True, "angles": 2, "conduction_only": True,
    },
]

# `iso_stiffness(k, μ)` takes *physical* moduli while the raw `TensISO{3}(a, b)`
# constructor takes `(3k, 2μ)`. Only the former is offered, which is how that
# trap is removed rather than merely documented.
PROPERTIES = [
    {
        "name": "iso_kmu", "label": "Isotropic (k, μ)", "order": 4,
        "builder": "iso_stiffness",
        "fields": [
            {"name": "k", "label": "k (bulk)", "type": "number", "default": 10.0},
            {"name": "mu", "label": "μ (shear)", "type": "number", "default": 5.0},
        ],
    },
    {
        "name": "iso_Enu", "label": "Isotropic (E, ν)", "order": 4,
        "builder": "iso_stiffness_E_nu",
        "fields": [
            {"name": "E", "label": "E (Young)", "type": "number", "default": 30.0},
            {"name": "nu", "label": "ν (Poisson)", "type": "number", "default": 0.2},
        ],
    },
    {
        "name": "ti_hoenig", "label": "Transversely isotropic (Hoenig)", "order": 4,
        "builder": "hoenig_stiffness",
        "fields": [
            {"name": "E1", "label": "E₁", "type": "number", "default": 30.0},
            {"name": "h", "label": "h", "type": "number", "default": 1.0},
            {"name": "nu1", "label": "ν₁", "type": "number", "default": 0.2},
            {"name": "nu2", "label": "ν₂", "type": "number", "default": 0.2},
            {"name": "gamma", "label": "γ", "type": "number", "default": 1.0},
        ],
    },
    {
        "name": "iso_conduction", "label": "Isotropic conductivity", "order": 2,
        "builder": "TensISO{2, 3}",
        "fields": [{"name": "k", "label": "κ", "type": "number", "default": 1.0}],
    },
    {
        "name": "void", "label": "Void / pore (near-zero)", "order": 4,
        "builder": "iso_stiffness",
        "doc": "Kept slightly non-zero so the tensor stays invertible.",
        "fields": [
            {"name": "k", "label": "k", "type": "number", "default": 1.0e-6},
            {"name": "mu", "label": "μ", "type": "number", "default": 1.0e-6},
        ],
    },
]

# `symmetrize` is an exact rotational average applied inside the kernel;
# `best_fit_*` is a least-squares projection used for reporting. Conflating
# them changes the numbers, so they are two separate lists.
SYMMETRIZE = [
    {"name": "none", "label": "None", "emit": "nothing"},
    {
        "name": "iso", "label": "Isotropic orientation average",
        "emit": "IsoSymmetrize()",
        "doc": "Uniform distribution of orientations, averaged exactly in the kernel.",
    },
    {
        "name": "ti", "label": "Transversely isotropic average",
        "emit": "TISymmetrize()",
        "doc": "Orientations uniformly distributed about an axis.",
    },
]

PROJECTIONS = [
    {"name": "none", "label": "As computed", "emit": None},
    {"name": "iso", "label": "Best isotropic fit", "emit": "best_fit_iso"},
    {"name": "ti", "label": "Best TI fit", "emit": "best_fit_ti"},
    {"name": "ortho", "label": "Best orthotropic fit", "emit": "best_fit_ortho"},
]

INTERFACES = [
    {"name": "PerfectInterface", "label": "Perfect", "fields": []},
    {
        "name": "SpringInterface", "label": "Spring (displacement jump)",
        "fields": [
            {"name": "kn", "label": "kₙ", "type": "number", "default": 1.0},
            {"name": "kt", "label": "kₜ", "type": "number", "default": 1.0},
        ],
    },
    {
        "name": "MembraneInterface", "label": "Membrane (traction jump)",
        "fields": [
            {"name": "k2D", "label": "k₂D", "type": "number", "default": 1.0},
            {"name": "mu2D", "label": "μ₂D", "type": "number", "default": 1.0},
        ],
    },
    {
        "name": "KapitzaInterface", "label": "Kapitza (thermal)",
        "fields": [{"name": "h", "label": "h", "type": "number", "default": 1.0}],
        "order": 2,
    },
    {
        "name": "SurfaceConductiveInterface", "label": "Surface conductive",
        "fields": [{"name": "ks", "label": "κₛ", "type": "number", "default": 1.0}],
        "order": 2,
    },
]

LENSES = [
    {
        "name": "amount", "label": "Phase amount (fraction or density)",
        "args": ["phase"],
        "doc": "The matrix amount is derived (1 − Σ f) and cannot be set.",
    },
    {"name": "property", "label": "Property component", "args": ["phase", "property", "index"]},
    {
        "name": "geometry", "label": "Geometry field",
        "args": ["phase", "field", "index"],
        "doc": "Changing a semi-axis reclassifies the shape trait.",
    },
    {"name": "shape_param", "label": "Distribution shape", "args": ["field", "index"]},
    {
        "name": "nested", "label": "Through a nested scale",
        "args": ["member", "property", "inner"],
        "doc": "Reaches into a Homogenized inner cell; crosses scales for AD.",
    },
]

VISCO = [
    {
        "name": "maxwell_iso", "label": "Maxwell (relaxation)", "mode": "relaxation",
        "fields": [
            {"name": "k", "label": "k", "type": "number", "default": 10.0},
            {"name": "mu", "label": "μ", "type": "number", "default": 5.0},
            {"name": "tau", "label": "τ (relaxation time)", "type": "number", "default": 1.0},
        ],
    },
    {
        "name": "kelvin_iso", "label": "Kelvin-Voigt (creep)", "mode": "creep",
        "fields": [
            {"name": "k", "label": "k", "type": "number", "default": 10.0},
            {"name": "mu", "label": "μ", "type": "number", "default": 5.0},
            {"name": "tau", "label": "τ (creep time)", "type": "number", "default": 1.0},
        ],
    },
    {
        "name": "heaviside", "label": "Elastic (Heaviside)", "mode": "relaxation",
        "fields": [
            {"name": "k", "label": "k", "type": "number", "default": 10.0},
            {"name": "mu", "label": "μ", "type": "number", "default": 5.0},
        ],
    },
    {
        "name": "custom", "label": "Custom J(t, t′) or R(t, t′)", "mode": "creep",
        "fields": [{
            "name": "expr", "label": "Julia expression in t, t′", "type": "code",
            "default": "1.0 / 10.0 * (1 + log(1 + (t - t′)))",
        }],
    },
]

#: Used until the sidecar answers. Replaced wholesale, never merged, so a
#: scheme MeanFieldHom drops does not linger in the interface.
FALLBACK_SCHEMES = [
    {"name": n, "options": [], "singleton": True}
    for n in (
        "AsymmetricSelfConsistent", "DifferentialScheme", "Dilute", "DiluteDual",
        "Laminated", "Maxwell", "MoriTanaka", "PonteCastanedaWillis", "Reuss",
        "SelfConsistent", "Voigt",
    )
]


def base_catalog() -> dict:
    """Everything the interface can offer without Julia."""
    return {
        "mfh_version": None,
        "julia_version": None,
        "introspected": False,
        "schemes": FALLBACK_SCHEMES,
        "geometries": GEOMETRIES,
        "properties": PROPERTIES,
        "symmetrize": SYMMETRIZE,
        "projections": PROJECTIONS,
        "interfaces": INTERFACES,
        "lenses": LENSES,
        "visco": VISCO,
        "hill_methods": ["auto", "analytical", "residues", "decuhr", "nestedquadgk"],
        "constraints": {
            "no_multiscale_in_alv": True,
            "matrix_amount_is_derived": True,
        },
    }


def merge(introspected: dict | None) -> dict:
    """The form definitions, with the live facts folded in when available."""
    cat = base_catalog()
    if not introspected:
        return cat
    cat["schemes"] = introspected.get("schemes") or cat["schemes"]
    cat["mfh_version"] = introspected.get("mfh_version")
    cat["julia_version"] = introspected.get("julia_version")
    cat["introspected"] = True
    if introspected.get("constraints"):
        cat["constraints"].update(introspected["constraints"])
    return cat
