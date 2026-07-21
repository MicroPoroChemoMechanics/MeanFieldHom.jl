# API — Schemes

Public types and functions of `MeanFieldHom.Schemes`.

## RVE / Phase / Amount

```@docs
RVE
Phase
AbstractAmount
VolumeFraction
CrackDensity
AbstractDistributionShape
UniformDistribution
AbstractSymmetrize
NoSymmetrize
IsoSymmetrize
TISymmetrize
phase_symmetrize
add_matrix!
add_phase!
matrix_phase
inclusion_phase_names
phase_property
matrix_property
volume_fraction
crack_density
matrix_volume_fraction
validate_rve
```

## Schemes

```@docs
HomogenizationScheme
Voigt
Reuss
Dilute
DiluteDual
MoriTanaka
Maxwell
PonteCastanedaWillis
SelfConsistent
AsymmetricSelfConsistent
DifferentialScheme
DifferentialTrajectory
Proportional
Sequential
CustomPath
Path
AndersonDefault
NewtonDefault
```

## Entry point

```@docs
homogenize
MeanFieldHom.Schemes.SCHEME_ALIAS
```
