# Changelog

## v0.2.0 — alignment with TensND 0.2 (breaking)

Follow-up to TensND 0.2's API unification. MeanFieldHom is iso-functional —
all outputs are unchanged — but every mention of a TensND symbol now uses
the new snake_case + UPPERCASE-acronym convention.

### Breaking changes

- `TensND.TensWalpole` references (type annotations, dispatch rules,
  constructor calls) now use `TensND.TensTI{4}`.  The struct layout is
  identical so numerical behaviour is unchanged.
- Accessor renames propagated from TensND: `getbasis` → `get_basis`,
  `tensbasis` → `tens_basis`, `invKM` → `inv_KM`, `getdata` → `get_data`,
  `getarray` → `get_array`, `getvar` → `get_var`, `getdim` → `get_dim`,
  `getorder` → `get_order`.
- Predicate renames: `isISO` → `is_ISO`, `isTI` → `is_TI`,
  `isOrtho` → `is_ORTHO`.
- Tensor factory renames in scripts and docs: `tensId2` → `tens_Id2`,
  `tensJ4` → `tens_J4`, `tensTI` → `tens_TI`, etc.

### Additions

None — functional surface unchanged.

### Migration guide

If you have your own code depending on MeanFieldHom dispatch, apply the
same renames as listed in TensND's v0.2 changelog. All MeanFieldHom tests
(2865) pass without behavioural change after migration.
