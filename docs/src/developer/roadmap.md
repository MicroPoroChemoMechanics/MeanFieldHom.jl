# Roadmap

- Multi-layer inclusions (concentric spheres, coated cylinders) via
  [`AbstractLayeredInclusion`](@ref).
- Mean-field schemes: dilute, Mori–Tanaka, self-consistent,
  Ponte-Castañeda–Willis, differential.
- Representative volume element (RVE) assembly and effective-property
  pipelines mirroring the C++ ECHOES `rve.h`.
- Viscoelastic constitutive laws (time-domain and Laplace–Carson).
- User-defined inclusions / algorithms via the open `_kernel` table.
- C++ <-> Julia bridge (CxxWrap / ccall) against the Blitz++-backed
  ECHOES core library.
