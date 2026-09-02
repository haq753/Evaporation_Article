# Physics verification

This document is the Gate A / Gate B deliverable: a record of every published worked
example the physics core has been checked against, with the comparison numbers, so a
reviewer (or a future session) does not have to re-derive them.

Two rules were followed throughout: every "expected" value below was independently
computed from the cited published formula or transcribed from a published worked
example, never taken from recollection; and a disagreement was diagnosed before being
treated as passing (see the Rn/Rnl notes below for the one case where a gap turned out
to be a legitimate published-convention difference rather than a bug).

## `ETo_FAO56()` — FAO-56 (Allen et al., 1998), Examples 18 and 20

| Quantity | Published | Computed | Note |
|---|---|---|---|
| Example 18 (Uccle, Brussels, 6 Jul), `ea` | 1.409 kPa | 1.409 | intermediate |
| Example 18, `Rs` | 22.07 MJ m⁻² d⁻¹ | 22.07 | intermediate |
| Example 18, `es` | 1.997 kPa | 1.997 | intermediate |
| Example 18, `VPD` | 0.589 kPa | 0.589 | intermediate |
| Example 18, `Rn` | 13.28 MJ m⁻² d⁻¹ | 13.28 | intermediate |
| Example 18, `ETo` | 3.9 mm/day | 3.8796 | final |
| Example 20 (Lyon, missing data), `ETo` | 4.6 mm/day | 4.5604 | final |

## Shared radiation geometry — `Ra_calc()`

| Site / date | Published or independently derived | Computed | Source |
|---|---|---|---|
| Uccle, 50.80°N, 6 Jul | 41.088 MJ m⁻² d⁻¹ | 41.088 | FAO-56 |
| Lyon, 45.72°N, 15 Jul | 40.555 MJ m⁻² d⁻¹ | 40.555 | FAO-56 |
| Equator, equinox | 37.824 MJ m⁻² d⁻¹ | 37.824 | Eq. 21-25, independently computed |
| Atatürk site bracket, 37.5°N, 21 Jun | 41.779 MJ m⁻² d⁻¹ | 41.779 | Eq. 21-25, independently computed |
| Atatürk site bracket, 37.5°N, 21 Dec | 15.069 MJ m⁻² d⁻¹ | 15.069 | Eq. 21-25, independently computed |

The equator and Atatürk-bracket rows replace figures that appeared in earlier project
notes (≈36.2 and ≈41.5 respectively) and that turned out to be wrong when actually
computed — not a defect in `Ra_calc()`, a defect in an unverified recollection. See
`PROJECT_HANDOVER.md` section 5.

## `Eopen_Penman()` — McMahon et al. (2013), *HESS* 17, 1331–1363, Supplement S19

**Worked Examples 1–3.** Station: Alice Springs Airport, Australia (23.7951°S, 546 m),
20 July 1980. Tmax 21.0°C, Tmin 2.0°C, RHmax 71%, RHmin 25%, wind run 51 km/day at 2 m
(u2 = 0.5903 m/s), sunshine 10.7 h (Rs = 17.1940 MJ m⁻² d⁻¹ as computed by the paper).
`G = 0`, `Tw_C = NULL` (no heat storage — this verifies the base combination equation,
not the heat-storage extension, which is Stage D).

| Quantity | Published (S19) | Computed | Agreement |
|---|---|---|---|
| Tmean | 11.5 °C | 11.500 | exact |
| es (mean of es(Tmax), es(Tmin)) | 1.5963 kPa | 1.59632 | ~1e-5 |
| ea | 0.5614 kPa | 0.56138 | ~1e-5 |
| VPD | 1.0349 kPa | 1.03495 | ~1e-5 |
| Ra | 23.6182 MJ m⁻² d⁻¹ | 23.61822 | ~1e-5 |
| Rso | 17.9716 MJ m⁻² d⁻¹ | 17.97158 | ~1e-5 |
| Rns (water, α=0.08) | 15.8184 MJ m⁻² d⁻¹ | 15.81848 | ~1e-5 |
| Rnl | 7.1784 MJ m⁻² d⁻¹ | 7.17437 | 0.06% |
| Rn (water) | 8.6401 MJ m⁻² d⁻¹ | 8.64411 | 0.05% |
| Atmospheric pressure | 95.01027 kPa | 95.01027 | exact |
| γ | 0.0632 kPa/°C | 0.06318 | ~1e-5 |
| Δ | 0.0898 kPa/°C | 0.08984 | ~1e-4 |
| EPenOW, McMahon's own wind function (`wind_fn = "penman1956"`) | 2.9797 mm/day | 2.9808 | 0.04% |

The Rnl/Rn gap (0.05–0.06%) is fully explained by a documented convention difference:
FAO-56 Eq. 39 uses `T + 273.16` (preserved exactly in `utils_physics.R`, per its header
note), while McMahon's Eq. S3.3 uses `T + 273.2`. Recomputing Rnl with 273.2 instead of
273.16 gives 7.17839, matching the published 7.1784 to five significant figures. Both
conventions are legitimate published choices; this is not a defect.

The final-answer gap (0.04%) is within the rounding cascade expected from a hand-worked
example that itself carries each intermediate to 4-5 significant figures — every
intermediate above already agrees to 1e-5, so the combination equation itself
(`rad_term`, `aero_term`, and their assembly) is confirmed correct.

### Wind function: two published formulas, one deliberate choice

`Eopen_Penman()`'s aerodynamic term can use either of two published wind functions,
selected with `wind_fn`:

- `"fao56_shuttleworth"` (**default**) — Shuttleworth's (1993) recasting of the
  land-calibrated FAO/Penman wind function, `6.43(1+0.536 u2)` MJ m⁻² d⁻¹ kPa⁻¹.
- `"penman1956"` — Penman's (1956) own open-water wind function,
  `f(u) = 1.313 + 1.381 u2` mm day⁻¹ kPa⁻¹, used throughout McMahon et al. (2013)'s
  worked examples (their Eq. S4.3).

For the Alice Springs case, these two published parameterisations give Eopen of 3.548
and 2.981 mm/day respectively — a **19% difference from wind-function choice alone**,
holding every other input identical. This is larger than the ~5% previously anticipated
for the separate question of wind-height/roughness conversion (10 m → 2 m; land vs.
open-water log-law), and needs its own line in Methods, not just a footnote.

**Decision (Aziz, 2 Sep 2026):** `fao56_shuttleworth` remains the default for the main
analysis — both formulas are published and neither is more "correct" a priori.
`wind_fn = "penman1956"` is kept as an explicit, one-argument option, used here to
reproduce the Gate B worked example and available as a documented sensitivity run.

## Test coverage

All of the above (plus the wind-function sensitivity as a pinned regression value) is
encoded in `tests/test_physics.R`: 17 `test_that()` blocks, 60 expectations, all passing
as of the commit that introduced this file. Re-run `testthat::test_dir("tests")` after
any change to `R/utils_physics.R` or `R/Penman_Combination_Open_ETo.R`.
