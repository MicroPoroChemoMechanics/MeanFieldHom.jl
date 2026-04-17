# Theoretical overview

`MeanFieldHom` implements the Fourier-space / Green-operator machinery
underpinning Eshelby-type inclusion problems and flat-crack analyses:

```math
P = \int_{\lVert\boldsymbol\xi\rVert=1} \Gamma(\boldsymbol\xi;\mathbb C_0)\,\mathrm d S
```

for ellipsoidal inclusions, and the line-integrated Green kernel
``\hat{\mathbf Q}^\star_{nn}`` for flat cracks.  Detailed derivations
are deferred to later versions of this manual; see the references
cited in the module docstrings and the companion paper.
