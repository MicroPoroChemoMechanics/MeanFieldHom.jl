# Viscoelastic composites

Every scheme seen so far takes an elastic stiffness tensor and returns
an elastic stiffness tensor. The **correspondence principle** extends
this, almost for free, to linear viscoelasticity: replace real moduli
by their complex, frequency-dependent counterparts, and the exact same
schemes carry through.

## Complex moduli in the frequency domain

A Maxwell-type viscoelastic shear modulus has the classical form

```math
G^*(\omega) = G_\infty + G_d\,\frac{i\omega\tau}{1+i\omega\tau},
```

real (elastic) at ``\omega \to 0`` and again at ``\omega \to \infty``,
complex in between. Building an RVE with complex-valued
`iso_stiffness` and calling [`homogenize`](@ref) exactly as
before propagates that complex modulus through the homogenization
algebra — every scheme is `Complex{Float64}`-safe:

```@example tutvisco
using MeanFieldHom
using TensND
using LinearAlgebra
using Plots
gr()  # headless backend; GKSwstype is set to "100" in make.jl

maxwell_G(ω; G_inf, G_d, τ = 1.0) = G_inf + G_d * (im * ω * τ) / (1 + im * ω * τ)

const f_inc = 0.3
ωs = exp10.(range(-2, 2; length = 40))

function C_eff_at(ω, scheme)
    G_m = maxwell_G(ω; G_inf = 10.0, G_d = 5.0)
    G_i = maxwell_G(ω; G_inf = 30.0, G_d = 15.0)
    r = RVE(:M; T = ComplexF64)
    add_matrix!(r, Ellipsoid(1.0), Dict(:C => iso_stiffness(ComplexF64(30.0), G_m)))
    add_phase!(r, :I, Ellipsoid(1.0), Dict(:C => iso_stiffness(ComplexF64(80.0), G_i)); fraction = ComplexF64(f_inc))
    _, μ = k_mu(homogenize(r, scheme, :C))
    return μ
end

y_mt = [C_eff_at(ω, MoriTanaka()) for ω in ωs]
y_sc = [C_eff_at(ω, SelfConsistent(; abstol = 1.0e-10, maxiters = 200)) for ω in ωs]

plt = plot(;
    xlabel = "ω", ylabel = "μ_eff", xscale = :log10,
    legend = :topleft, framestyle = :box, size = (760, 480),
)
plot!(plt, ωs, real.(y_mt); label = "storage — MT", color = :blue, lw = 2)
plot!(plt, ωs, imag.(y_mt); label = "loss — MT", color = :blue, lw = 2, ls = :dash)
plot!(plt, ωs, real.(y_sc); label = "storage — SC", color = :red, lw = 2)
plot!(plt, ωs, imag.(y_sc); label = "loss — SC", color = :red, lw = 2, ls = :dash)
plt
```

The real part (storage modulus) rises from the low-frequency to the
high-frequency elastic limit, and the imaginary part (loss modulus)
peaks in between, at the frequency where the microstructure's relaxation
time and the loading period coincide — the qualitative signature of any
viscoelastic composite, MT and SC differing only in how strongly the
inclusion phase influences the peak's height and position through
interaction.

!!! note "Note the RVE's element type"
    Complex-valued homogenization needs `RVE(:M; T = ComplexF64)`; the
    default `RVE(:M)` fixes the volume-fraction element type to
    `Float64` and does not accept complex phase properties.

## A first taste of time-domain ageing viscoelasticity

The frequency-domain view above assumes properties that do not evolve
with the material's age. For **ageing** viscoelasticity — a material
whose relaxation spectrum itself changes with time, as in curing cement
paste — `MeanFieldHom` provides a full time-domain (ALV) pipeline built
on discretized Volterra operators. A [`ViscoLaw`](@ref) wraps a
two-argument relaxation kernel `R(t, t')`, and
[`homogenize_alv`](@ref) takes the place of `homogenize`:

```@example tutvisco
function R_iso(t, tp)
    α = 3 * (3.0 + 2.0 * exp(-(t - tp) / 1.0))
    β = 2 * (1.0 + 1.0 * exp(-(t - tp) / 0.5))
    return TensISO{3}(α, β)
end
law_M = ViscoLaw(R_iso, :relaxation)

rve = RVE(:M)
add_matrix!(rve, Ellipsoid(1.0), Dict(:C => law_M))
add_phase!(rve, :I, Ellipsoid(1.0), Dict(:C => heaviside_law(iso_stiffness(60.0, 20.0))); fraction = 0.2)

times = collect(range(0.0, 5.0; length = 20))
C_eff = homogenize_alv(rve, MoriTanaka(), :C; times = times)
size(C_eff)
```

The result is a discretized Volterra operator (a dense block matrix,
here `40×40` since each of the 20 time steps carries a `2×2` iso
block), not a single tensor — reading effective moduli back out of it,
handling cracks in ALV, and differentiating through the whole pipeline
are covered in full in the
[Viscoelasticity manual](../manual/viscoelasticity.md). The next two
tutorials return to elastic problems and build up to differentiating
`homogenize` itself.
