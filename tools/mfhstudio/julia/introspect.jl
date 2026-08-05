# =============================================================================
#  introspect.jl — the MFH feature catalogue, discovered at run time.
#
#  Hard-coding the list of schemes, interfaces and inclusion types in the web
#  UI would make the interface silently fall behind MeanFieldHom. Everything
#  the UI offers is therefore enumerated from the loaded package, in the same
#  spirit as `echoes2mfh check-drift`.
# =============================================================================

using MeanFieldHom
using TensND
using InteractiveUtils: subtypes

"""
    concrete_subtypes(T) -> Vector{Type}

Every instantiable leaf under `T`, depth-first. Abstract intermediate layers
are traversed but not reported.
"""
function concrete_subtypes(T::Type)
    out = Type[]
    for S in subtypes(T)
        if isabstracttype(S)
            append!(out, concrete_subtypes(S))
        else
            push!(out, S)
        end
    end
    return out
end

# ── Constructor keywords ─────────────────────────────────────────────────────
#
# The scheme constructors carry the solver options (`abstol`, `maxiters`,
# `select_best`, `nsteps`, …). Reading them off the methods keeps the UI's
# option list exact rather than approximate.

function ctor_keywords(T::Type)
    names = Set{Symbol}()
    for m in methods(T)
        for kw in Base.kwarg_decl(m)
            s = String(kw)
            # `kwargs...` shows up as a trailing `...` name
            endswith(s, "...") && continue
            push!(names, kw)
        end
    end
    return sort!(collect(names); by = String)
end

# Several schemes take their solver options through a `kwargs...` bag
# (`SelfConsistent(; algorithm = …, kwargs...)`), so `kwarg_decl` cannot see
# them — and the bag accepts *anything*, so probing the constructor proves
# nothing either: `SelfConsistent(; nsteps = 3)` succeeds while `nsteps` is
# meaningless there.
#
# What the schemes do publish is the list of keys they actually read. Reading
# those constants keeps the interface exactly in step with the schemes, and
# the fallback below means a rename degrades to "no options offered" rather
# than to a wrong list.
const OPTION_KEY_CONSTANTS = Dict(
    "SelfConsistent" => :_SC_SOLVER_KWARGS,
    "AsymmetricSelfConsistent" => :_SC_SOLVER_KWARGS,
    "DifferentialScheme" => :_DIFF_RESERVED_OPTIONS,
)

const OPTION_DEFAULTS = Dict(
    :abstol => 1.0e-8, :reltol => 1.0e-6, :maxiters => 100,
    :damping => 0.0, :select_best => false, :verbose => false,
    :nsteps => 100, :formulation => "stiffness",
)

"""
    consumed_options(S) -> Vector{Symbol}

The option keys scheme `S` actually reads, taken from the constant it
declares for the purpose. Empty when the scheme declares none.
"""
function consumed_options(S::Type)
    name = String(nameof(S))
    key = get(OPTION_KEY_CONSTANTS, name, nothing)
    key === nothing && return Symbol[]
    mod = parentmodule(S)
    isdefined(mod, key) || return Symbol[]
    return collect(Symbol, getfield(mod, key))
end

"""
Options that hold an object rather than a number (`algorithm`, `trajectory`,
`alg`). They are reported so the UI can display the default, but not offered
as free-form inputs.
"""
const OPAQUE_OPTIONS = Set([:algorithm, :trajectory, :alg])

function scheme_entry(S::Type)
    consumed = consumed_options(S)
    declared = ctor_keywords(S)

    opts = Dict{String, Any}[]
    for k in consumed
        push!(
            opts, Dict(
                "name" => String(k),
                "default" => get(OPTION_DEFAULTS, k, nothing),
                "editable" => !(k in OPAQUE_OPTIONS),
            )
        )
    end
    for k in declared
        k in consumed && continue
        push!(
            opts, Dict(
                "name" => String(k), "default" => nothing,
                "editable" => !(k in OPAQUE_OPTIONS),
            )
        )
    end
    sort!(opts; by = d -> d["name"])

    return Dict(
        "name" => String(nameof(S)),
        "options" => opts,
        # A scheme with nothing to configure is a singleton: `Voigt()`.
        "singleton" => all(d -> !d["editable"], opts),
    )
end

# ── Geometry ────────────────────────────────────────────────────────────────
#
# Each entry describes what the UI must ask for. `fields` drives the form; the
# emitted Julia is built by the Python side from the same names.

const GEOMETRY_FORMS = [
    Dict(
        "name" => "Spheroid", "kind" => "spheroid", "dim" => 3,
        "doc" => "Axisymmetric ellipsoid; ω < 1 oblate, ω > 1 prolate, ω = 1 sphere.",
        "fields" => [
            Dict("name" => "omega", "label" => "ω (aspect ratio)", "type" => "number", "default" => 1.0),
        ],
        "angles" => 2,
    ),
    Dict(
        "name" => "Ellipsoid", "kind" => "ellipsoid", "dim" => 3,
        "doc" => "General ellipsoid of semi-axes (a, b, c).",
        "fields" => [
            Dict("name" => "a", "label" => "a", "type" => "number", "default" => 1.0),
            Dict("name" => "b", "label" => "b", "type" => "number", "default" => 1.0),
            Dict("name" => "c", "label" => "c", "type" => "number", "default" => 1.0),
        ],
        "angles" => 3,
    ),
    Dict(
        "name" => "Cylinder", "kind" => "cylinder", "dim" => 3,
        "doc" => "Infinite elliptic cylinder of transverse semi-axes (b, c).",
        "fields" => [
            Dict("name" => "b", "label" => "b", "type" => "number", "default" => 1.0),
            Dict("name" => "c", "label" => "c", "type" => "number", "default" => 1.0),
        ],
        "angles" => 3,
    ),
    Dict(
        "name" => "PennyCrack", "kind" => "penny_crack", "dim" => 3,
        "doc" => "Circular flat crack of radius a. Enters the RVE with a crack density.",
        "fields" => [
            Dict("name" => "a", "label" => "a (radius)", "type" => "number", "default" => 1.0),
        ],
        "angles" => 2, "amount" => "density",
    ),
    Dict(
        "name" => "EllipticCrack", "kind" => "elliptic_crack", "dim" => 3,
        "doc" => "Elliptic flat crack, semi-axes a ≥ b.",
        "fields" => [
            Dict("name" => "a", "label" => "a", "type" => "number", "default" => 1.0),
            Dict("name" => "b", "label" => "b", "type" => "number", "default" => 0.5),
        ],
        "angles" => 3, "amount" => "density",
    ),
    Dict(
        "name" => "RibbonCrack", "kind" => "ribbon_crack", "dim" => 3,
        "doc" => "Tunnel crack of half-width b, unbounded along its length.",
        "fields" => [
            Dict("name" => "b", "label" => "b (half-width)", "type" => "number", "default" => 1.0),
        ],
        "angles" => 3, "amount" => "density",
    ),
    Dict(
        "name" => "LayeredSphere", "kind" => "layered_sphere", "dim" => 3,
        "doc" => "Concentric layers, ascending radii, r = 0 implicit at the centre.",
        "fields" => [], "layered" => true, "angles" => 0,
    ),
    Dict(
        "name" => "LayeredSpheroid", "kind" => "layered_spheroid", "dim" => 3,
        "doc" => "Confocal spheroidal layers. Conduction only.",
        "fields" => [
            Dict("name" => "Nseries", "label" => "N (series order)", "type" => "integer", "default" => 5),
        ],
        "layered" => true, "angles" => 2, "conduction_only" => true,
    ),
]

# ── Material properties ─────────────────────────────────────────────────────
#
# `iso_stiffness(k, μ)` takes *physical* moduli while the raw `TensISO{3}(a, b)`
# constructor takes `(3k, 2μ)`. The UI only ever offers the former, which is
# how the trap is removed rather than merely documented.

const PROPERTY_FORMS = [
    Dict(
        "name" => "iso_kmu", "label" => "Isotropic (k, μ)", "order" => 4,
        "builder" => "iso_stiffness",
        "fields" => [
            Dict("name" => "k", "label" => "k (bulk)", "type" => "number", "default" => 10.0),
            Dict("name" => "mu", "label" => "μ (shear)", "type" => "number", "default" => 5.0),
        ],
    ),
    Dict(
        "name" => "iso_Enu", "label" => "Isotropic (E, ν)", "order" => 4,
        "builder" => "iso_stiffness_E_nu",
        "fields" => [
            Dict("name" => "E", "label" => "E (Young)", "type" => "number", "default" => 30.0),
            Dict("name" => "nu", "label" => "ν (Poisson)", "type" => "number", "default" => 0.2),
        ],
    ),
    Dict(
        "name" => "ti_hoenig", "label" => "Transversely isotropic (Hoenig)", "order" => 4,
        "builder" => "hoenig_stiffness",
        "fields" => [
            Dict("name" => "E1", "label" => "E₁", "type" => "number", "default" => 30.0),
            Dict("name" => "h", "label" => "h", "type" => "number", "default" => 1.0),
            Dict("name" => "nu1", "label" => "ν₁", "type" => "number", "default" => 0.2),
            Dict("name" => "nu2", "label" => "ν₂", "type" => "number", "default" => 0.2),
            Dict("name" => "gamma", "label" => "γ", "type" => "number", "default" => 1.0),
        ],
        "axis" => true,
    ),
    Dict(
        "name" => "iso_conduction", "label" => "Isotropic conductivity", "order" => 2,
        "builder" => "TensISO{2, 3}",
        "fields" => [
            Dict("name" => "k", "label" => "κ", "type" => "number", "default" => 1.0),
        ],
    ),
    Dict(
        "name" => "void", "label" => "Void / pore (near-zero)", "order" => 4,
        "builder" => "iso_stiffness",
        "doc" => "Kept slightly non-zero so the tensor stays invertible.",
        "fields" => [
            Dict("name" => "k", "label" => "k", "type" => "number", "default" => 1.0e-6),
            Dict("name" => "mu", "label" => "μ", "type" => "number", "default" => 1.0e-6),
        ],
    ),
]

# ── Symmetrization ──────────────────────────────────────────────────────────
#
# `symmetrize` is an exact rotational average applied inside the kernel;
# `best_fit_*` is a least-squares projection used for reporting. Conflating
# them changes the numbers, so the catalogue keeps them apart and says so.

const SYMMETRIZE_FORMS = [
    Dict("name" => "none", "label" => "None", "emit" => "nothing"),
    Dict(
        "name" => "iso", "label" => "Isotropic orientation average",
        "emit" => "IsoSymmetrize()",
        "doc" => "Uniform distribution of orientations, averaged exactly in the kernel.",
    ),
    Dict(
        "name" => "ti", "label" => "Transversely isotropic average",
        "emit" => "TISymmetrize()",
        "doc" => "Orientations uniformly distributed about an axis.",
    ),
]

const REPORT_PROJECTIONS = [
    Dict("name" => "none", "label" => "As computed", "emit" => nothing),
    Dict("name" => "iso", "label" => "Best isotropic fit", "emit" => "best_fit_iso"),
    Dict("name" => "ti", "label" => "Best TI fit", "emit" => "best_fit_ti"),
    Dict("name" => "ortho", "label" => "Best orthotropic fit", "emit" => "best_fit_ortho"),
]

# ── Interfaces (layered inclusions) ─────────────────────────────────────────

const INTERFACE_FORMS = [
    Dict("name" => "PerfectInterface", "label" => "Perfect", "fields" => []),
    Dict(
        "name" => "SpringInterface", "label" => "Spring (displacement jump)",
        "fields" => [
            Dict("name" => "kn", "label" => "kₙ", "type" => "number", "default" => 1.0),
            Dict("name" => "kt", "label" => "kₜ", "type" => "number", "default" => 1.0),
        ],
    ),
    Dict(
        "name" => "MembraneInterface", "label" => "Membrane (traction jump)",
        "fields" => [
            Dict("name" => "k2D", "label" => "k₂D", "type" => "number", "default" => 1.0),
            Dict("name" => "mu2D", "label" => "μ₂D", "type" => "number", "default" => 1.0),
        ],
    ),
    Dict(
        "name" => "KapitzaInterface", "label" => "Kapitza (thermal)",
        "fields" => [Dict("name" => "h", "label" => "h", "type" => "number", "default" => 1.0)],
        "order" => 2,
    ),
    Dict(
        "name" => "SurfaceConductiveInterface", "label" => "Surface conductive",
        "fields" => [Dict("name" => "ks", "label" => "κₛ", "type" => "number", "default" => 1.0)],
        "order" => 2,
    ),
]

# ── Sensitivity lenses ──────────────────────────────────────────────────────

const LENS_FORMS = [
    Dict(
        "name" => "amount", "label" => "Phase amount (fraction or density)",
        "args" => ["phase"],
        "doc" => "The matrix amount is derived (1 − Σ f) and cannot be set.",
    ),
    Dict(
        "name" => "property", "label" => "Property component",
        "args" => ["phase", "property", "index"],
    ),
    Dict(
        "name" => "geometry", "label" => "Geometry field",
        "args" => ["phase", "field", "index"],
        "doc" => "Changing a semi-axis reclassifies the shape trait.",
    ),
    Dict("name" => "shape_param", "label" => "Distribution shape", "args" => ["field", "index"]),
    Dict(
        "name" => "nested", "label" => "Through a nested scale",
        "args" => ["member", "property", "inner"],
        "doc" => "Reaches into a Homogenized inner cell; crosses scales for AD.",
    ),
]

# ── Viscoelasticity ─────────────────────────────────────────────────────────

const VISCO_FORMS = [
    Dict(
        "name" => "maxwell_iso", "label" => "Maxwell (relaxation)",
        "mode" => "relaxation",
        "fields" => [
            Dict("name" => "k", "label" => "k", "type" => "number", "default" => 10.0),
            Dict("name" => "mu", "label" => "μ", "type" => "number", "default" => 5.0),
            Dict("name" => "tau", "label" => "τ (relaxation time)", "type" => "number", "default" => 1.0),
        ],
    ),
    Dict(
        "name" => "kelvin_iso", "label" => "Kelvin-Voigt (creep)",
        "mode" => "creep",
        "fields" => [
            Dict("name" => "k", "label" => "k", "type" => "number", "default" => 10.0),
            Dict("name" => "mu", "label" => "μ", "type" => "number", "default" => 5.0),
            Dict("name" => "tau", "label" => "τ (creep time)", "type" => "number", "default" => 1.0),
        ],
    ),
    Dict(
        "name" => "heaviside", "label" => "Elastic (Heaviside)",
        "mode" => "relaxation",
        "fields" => [
            Dict("name" => "k", "label" => "k", "type" => "number", "default" => 10.0),
            Dict("name" => "mu", "label" => "μ", "type" => "number", "default" => 5.0),
        ],
    ),
    Dict(
        "name" => "custom", "label" => "Custom J(t, t′) or R(t, t′)",
        "mode" => "creep",
        "fields" => [
            Dict(
                "name" => "expr", "label" => "Julia expression in t, t′",
                "type" => "code", "default" => "1.0 / 10.0 * (1 + log(1 + (t - t′)))",
            ),
        ],
    ),
]

"""
    catalogue() -> Dict

The whole feature list the web UI is built from.
"""
function catalogue()
    schemes = [scheme_entry(S) for S in concrete_subtypes(MeanFieldHom.HomogenizationScheme)]
    sort!(schemes; by = d -> d["name"])
    return Dict(
        "mfh_version" => string(pkgversion(MeanFieldHom)),
        "julia_version" => string(VERSION),
        "schemes" => schemes,
        "geometries" => GEOMETRY_FORMS,
        "properties" => PROPERTY_FORMS,
        "symmetrize" => SYMMETRIZE_FORMS,
        "projections" => REPORT_PROJECTIONS,
        "interfaces" => INTERFACE_FORMS,
        "lenses" => LENS_FORMS,
        "visco" => VISCO_FORMS,
        "hill_methods" => ["auto", "analytical", "residues", "decuhr", "nestedquadgk"],
        # Documented constraint: an inner Homogenized cannot sit inside an
        # ageing-viscoelastic chain, because the inner result would have to be
        # re-expressible as a ViscoLaw.
        "constraints" => Dict(
            "no_multiscale_in_alv" => true,
            "matrix_amount_is_derived" => true,
        ),
    )
end
