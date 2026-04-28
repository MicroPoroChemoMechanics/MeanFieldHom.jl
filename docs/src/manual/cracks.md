# Cracks

```julia
using MeanFieldHom, TensND
E, ν = 210.0, 0.3
k = E/(3*(1-2ν)); μ = E/(2*(1+ν))
C₀ = TensISO{3}(3k, 2μ)

# Penny-shaped crack — size-independent COD tensor B
pc = PennyCrack(1.0)
B  = cod_tensor(pc, C₀)

# Size-independent compliance contribution tensor H = (3/4) n̂ ⊗ˢ B ⊗ˢ n̂
H  = compliance_contribution(pc, C₀)

# Dilute compliance correction ΔS from the Budiansky density ε³ᵈ = N a b²
ε³ᵈ = 0.05
ΔS  = delta_compliance(pc, H, ε³ᵈ)      # = (4π/3) ε³ᵈ H

# Ribbon crack — same pattern, ε²ᵈ = N b² and ΔS = π ε²ᵈ H
r   = RibbonCrack(0.5)
H_r = compliance_contribution(r, C₀)    # H = (2/π) n̂ ⊗ˢ B ⊗ˢ n̂
ΔS_r = delta_compliance(r, H_r, 0.05)

# Thermal / conductivity — scalar COD b and rank-1 resistivity tensor R
K₀ = TensISO{3}(1.0)
b  = cod_tensor(pc, K₀)                  # scalar
R  = compliance_contribution(pc, K₀)     # R = (3/4) b (ŵ⊗ŵ)
ΔR = delta_resistivity(pc, R, 0.05)      # = (4π/3) ε³ᵈ R
```

## Cracks with finite interface stiffness (Sevostianov)

A flat crack carrying a **spring-like interface elasticity** with
stiffness 2-tensor ``\\mathbf K`` (3 × 3 symmetric, e.g. iso with
``K_n`` normal + ``K_t`` tangential) modifies the COD via

```math
\\mathbf B_{\\text{eff}} = (b\\mathbf K + \\mathbf B^{-1})^{-1}
                          = \\mathbf B \\cdot (\\mathbf I + b\\mathbf K\\mathbf B)^{-1}
```

with `b` = `semi_minor(crack)`. Limits : ``\\mathbf K = 0`` →
traction-free (recovers ``\\mathbf B``); ``\\mathbf K \\to \\infty`` →
rigid bond (``\\mathbf B_{\\text{eff}} \\to 0``).

```julia
# Elasticity : iso interface stiffness K = 5·𝟏
B_eff = cod_tensor(pc, C₀; K_interface = TensISO{3}(5.0))
H_eff = compliance_contribution(pc, C₀; K_interface = TensISO{3}(5.0))

# Conductivity (Kapitza scalar interface conductance α)
b_eff = cod_tensor(pc, K₀; α_interface = 1.0)
R_eff = compliance_contribution(pc, K₀; α_interface = 1.0)
```

When building an `RVE` for a `homogenize` call, attach the interface
data as **phase properties** so the dispatcher picks them up
automatically :

```julia
rve = RVE(:M)
add_matrix!(rve, Ellipsoid(1.0), Dict(:C => C₀, :K => K₀))
add_phase!(rve, :CRACK, PennyCrack(1.0),
            Dict(:C => C₀,
                  :K_interface => TensISO{3}(5.0),  # elastic interface
                  :K => K₀,
                  :α_interface => 1.0);             # Kapitza scalar
            density = 0.10, symmetrize = :iso)

C_eff = homogenize(rve, MoriTanaka(), :C)
K_eff = homogenize(rve, MoriTanaka(), :K)
```

For SC-type schemes on cracked RVEs use
[`AsymmetricSelfConsistent`](@ref) — the symmetric `SelfConsistent`
form does not handle cracks (its strain-concentration tensor is
singular).

For the **time-dependent** (ALV) version with `Rn(t,t')` and
`Rt(t,t')` ageing interface kernels, see the
[Viscoelasticity manual](viscoelasticity.md#5-cracks-in-alv).
References : [@sevostianov2002], [@barthelemyIJES2019].
