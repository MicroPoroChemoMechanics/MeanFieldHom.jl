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
