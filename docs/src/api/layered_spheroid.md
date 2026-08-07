# [API — LayeredSpheroid](@id api-layered-spheroid)

`layer_count`, `layer_modulus`, `layer_interface` and
`layer_volume_fraction` are shared generics extended from
`LayeredSpheres` — see [API — LayeredSphere](layered_sphere.md)
for their docstrings; they apply unchanged to `LayeredSpheroid`.

```@docs
MeanFieldHomogenization.LayeredSpheroids
MeanFieldHomogenization.LayeredSpheroids.LayeredSpheroid
MeanFieldHomogenization.LayeredSpheroids.layered_spheroid_from_fractions
MeanFieldHomogenization.LayeredSpheroids.layer_q
MeanFieldHomogenization.LayeredSpheroids.layer_semiaxes
MeanFieldHomogenization.LayeredSpheroids.outer_semiaxes
MeanFieldHomogenization.LayeredSpheroids.spheroid_state_sequence
MeanFieldHomogenization.LayeredSpheroids.spheroid_ba_ratios
MeanFieldHomogenization.LayeredSpheroids.get_layer
MeanFieldHomogenization.LayeredSpheroids.local_temperature
MeanFieldHomogenization.LayeredSpheroids.local_gradient
MeanFieldHomogenization.LayeredSpheroids.local_flux
MeanFieldHomogenization.LayeredSpheroids.coupling_matrices
MeanFieldHomogenization.LayeredSpheroids.legendre_odd
```
