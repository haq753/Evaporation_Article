# ===============================================================
# tests/test_physics.R
#
# Run with:  testthat::test_file("tests/test_physics.R")
# Or from the project root:  testthat::test_dir("tests")
#
# Everything here is deterministic and needs no downloaded data.
# It should pass before any ERA5-Land file is opened, and again after
# any edit to R/utils_physics.R.
# ===============================================================

library(testthat)

source(here::here("R", "utils_physics.R"))
source(here::here("R", "FAO56_ETo.R"))
source(here::here("R", "Penman_Combination_Open_ETo.R"))
source(here::here("R", "Wind_10m_to_2m.R"))


# ---------------------------------------------------------------
# 1. Unit checks on the primitives
# ---------------------------------------------------------------

test_that("saturation vapour pressure matches published values", {
  expect_equal(es_kPa(31.0), 4.492590,  tolerance = 1e-5)
  expect_equal(es_kPa(20.0), 2.338, tolerance = 1e-3)
  expect_equal(es_kPa(0.0),  0.611, tolerance = 1e-3)
})

test_that("slope and psychrometric constant match FAO-56 tables", {
  expect_equal(slope_es(20.0),    0.1448, tolerance = 1e-3)
  expect_equal(gamma_psy(101.3),  0.0674, tolerance = 1e-3)
})

test_that("dewpoint round-trips through vapour pressure", {
  td <- 8.2
  expect_equal(td_from_ea(es_kPa(td)), td, tolerance = 1e-8)
  expect_equal(es_kPa(8.2), 1.087465, tolerance = 1e-5)
})

test_that("relative humidity is consistent with the Tetens es", {
  # T = 31 C, RH = 24 % -> ea ~ 1.08 kPa, VPD ~ 3.42 kPa
  ea  <- ea_from_rh(31.0, 24.0)
  expect_equal(ea, 1.080, tolerance = 5e-3)
  expect_equal(es_kPa(31.0) - ea, 3.42, tolerance = 5e-3)
  expect_equal(rh_from_td(31.0, td_from_ea(ea)), 24.0, tolerance = 1e-8)
})

test_that("extraterrestrial radiation matches FAO-56 Annex 2", {
  expect_equal(Ra_calc(50.80, 187), 41.088, tolerance = 1e-3)  # Uccle, 6 Jul
  expect_equal(Ra_calc(45.72, 196), 40.555, tolerance = 1e-3)  # Lyon, 15 Jul
  expect_equal(Ra_calc(0, 80),      37.82427,   tolerance = 1e-5)  # equator, equinox
})

test_that("extraterrestrial radiation is correct at the Ataturk site bracket", {
  # 37.5 N, summer and winter solstice. Independently computed from the
  # FAO-56 Eq. 21-25 Ra formula (not recalled): 41.7790 and 15.0687.
  # Earlier recalled figures of ~41.5 / ~16.1 were wrong; see handover
  # doc section 5. These are regression values pinned to the same
  # formula already verified against Uccle/Lyon/equator above.
  expect_equal(Ra_calc(37.5, 172), 41.779, tolerance = 1e-3)  # 21 Jun
  expect_equal(Ra_calc(37.5, 355), 15.069, tolerance = 1e-3)  # 21 Dec
})

test_that("daylight hours are sane and clamped", {
  expect_equal(daylight_hours(50.80, 187), 16.14, tolerance = 1e-2)
  expect_equal(daylight_hours(0, 80),      12.0,  tolerance = 1e-2)
  # Poleward of the Arctic circle in midsummer the clamp bites and
  # returns 24 h rather than NaN.
  expect_false(is.na(daylight_hours(80, 172)))
  expect_equal(daylight_hours(80, 172), 24.0, tolerance = 1e-6)
})

test_that("wind conversion matches both branches", {
  expect_equal(wind10_to_wind2(2.78),                          2.0793, tolerance = 1e-3)
  expect_equal(wind10_to_wind2(2.78, method = "log", z0 = 2e-4), 2.3665, tolerance = 1e-3)
})


# ---------------------------------------------------------------
# 2. FAO-56 Example 18 — Uccle (Brussels), 6 July
# ---------------------------------------------------------------
# Published inputs:  lat 50 deg 48' N, elev 100 m
#                    Tmax 21.5 C, Tmin 12.3 C
#                    RHmax 84 %, RHmin 63 %
#                    u10 = 10 km/h = 2.78 m/s
#                    actual sunshine n = 9.25 h
# Published answer:  ETo = 3.9 mm/day
#
# Conversions needed because the function takes Td and Rs directly:
#   ea  from RHmax/RHmin via FAO-56 Eq. 17, then inverted to Td
#   Rs  from n/N via the Angstrom relation, FAO-56 Eq. 35

test_that("FAO-56 Example 18 reproduces", {
  
  lat  <- 50.80
  elev <- 100
  J    <- 187
  Tmax <- 21.5
  Tmin <- 12.3
  
  # ea from RHmax and RHmin (FAO-56 Eq. 17)
  ea <- (es_kPa(Tmin) * 84 / 100 + es_kPa(Tmax) * 63 / 100) / 2
  expect_equal(ea, 1.409, tolerance = 1e-3)   # published intermediate
  
  # Rs from sunshine hours (FAO-56 Eq. 35, as = 0.25, bs = 0.50)
  Ra <- Ra_calc(lat, J)
  N  <- daylight_hours(lat, J)
  Rs <- (0.25 + 0.50 * 9.25 / N) * Ra
  expect_equal(Rs, 22.07, tolerance = 5e-3)   # published intermediate
  
  df <- data.frame(
    date           = as.Date("1998-07-06"),
    Tmax_C         = Tmax,
    Tmin_C         = Tmin,
    Tdmean_C       = td_from_ea(ea),
    u2             = wind10_to_wind2(2.78),
    ssrd_MJ_m2_day = Rs
  )
  
  out <- ETo_FAO56(df, lat_deg = lat, elev_m = elev, verbose = FALSE)
  
  # Published intermediates
  expect_equal(out$es,     1.997,  tolerance = 1e-3)
  expect_equal(out$vpd,    0.589,  tolerance = 5e-3)
  expect_equal(out$Rn_MJ, 13.28,   tolerance = 5e-3)
  
  # Published answer, given to 2 significant figures as 3.9
  expect_equal(out$ETo_mm, 3.88, tolerance = 5e-3)
  expect_equal(round(out$ETo_mm, 1), 3.9)
})


# ---------------------------------------------------------------
# 3. FAO-56 Example 20 — Lyon, July, with missing data
# ---------------------------------------------------------------
# Published inputs:  lat 45 deg 43' N, elev 200 m
#                    Tmax 26.6 C, Tmin 14.8 C (July monthly means)
#                    no humidity, radiation or wind measurements
# Assumptions in the published solution:
#                    Tdew = Tmin;  u2 = 2 m/s
#                    Rs = kRs * sqrt(Tmax - Tmin) * Ra, kRs = 0.16
# Published answer:  ETo = 4.6 mm/day

test_that("FAO-56 Example 20 reproduces", {
  
  lat  <- 45.72
  elev <- 200
  J    <- 196
  Tmax <- 26.6
  Tmin <- 14.8
  
  Ra <- Ra_calc(lat, J)
  Rs <- 0.16 * sqrt(Tmax - Tmin) * Ra
  expect_equal(Rs, 22.29, tolerance = 5e-3)   # published intermediate
  
  df <- data.frame(
    date           = as.Date("1998-07-15"),
    Tmax_C         = Tmax,
    Tmin_C         = Tmin,
    Tdmean_C       = Tmin,      # the missing-humidity assumption
    u2             = 2.0,
    ssrd_MJ_m2_day = Rs
  )
  
  out <- ETo_FAO56(df, lat_deg = lat, elev_m = elev, verbose = FALSE)
  
  expect_equal(out$es,   2.583, tolerance = 1e-3)
  expect_equal(out$ea,   1.68,  tolerance = 5e-3)
  expect_equal(out$vpd,  0.90,  tolerance = 1e-2)
  
  expect_equal(out$ETo_mm, 4.56, tolerance = 5e-3)
  expect_equal(round(out$ETo_mm, 1), 4.6)
})


# ---------------------------------------------------------------
# 4. Eopen_Penman Gate B — McMahon et al. (2013), HESS 17, 1331-1363,
#    Supplement Section S19, Worked Examples 1-3. Alice Springs
#    Airport, 20 Jul 1980, G = 0 (no heat storage).
# ---------------------------------------------------------------
# Published inputs: lat 23.7951 S, elev 546 m
#                    Tmax 21.0 C, Tmin 2.0 C, RHmax 71 %, RHmin 25 %
#                    wind run 51 km/day at 2 m -> u2 = 0.5903 m/s
#                    Rs (given) = 17.1940 MJ m-2 d-1
# Published answer:  EPenOW = 2.9797 mm/day, using the Penman (1956)
#                    wind function f(u) = 1.313 + 1.381 u2 (Eq. S4.3) -
#                    NOT this project's default wind function, so the
#                    test below calls wind_fn = "penman1956" explicitly.
#
# ea is supplied here as an equivalent dewpoint (prepare_forcing() takes
# Tdmean_C, not RH), inverted from the paper's own RHmax/RHmin route
# (Eq. S2.7) so the input is equivalent, not independent, for that one
# step.
#
# Rnl/Rn are expected to differ from the paper by ~0.06%: FAO-56 Eq. 39
# uses T + 273.16, McMahon uses T + 273.2 (see utils_physics.R header).
# Both conventions are legitimate; the difference is not a bug.

test_that("Eopen_Penman reproduces McMahon et al. (2013) Worked Example 3", {

  lat  <- -23.7951
  elev <- 546
  Tmax <- 21.0
  Tmin <- 2.0
  u2   <- 51 * 1000 / (24 * 60 * 60)   # wind run -> m/s
  expect_equal(u2, 0.5903, tolerance = 1e-3)   # published intermediate

  # ea from RHmax/RHmin (Eq. S2.7), same form as FAO-56 Eq. 17
  ea <- (es_kPa(Tmin) * 71 / 100 + es_kPa(Tmax) * 25 / 100) / 2
  expect_equal(ea, 0.5614, tolerance = 1e-3)   # published intermediate

  df <- data.frame(
    date           = as.Date("1980-07-20"),
    Tmax_C         = Tmax,
    Tmin_C         = Tmin,
    Tdmean_C       = td_from_ea(ea),
    u2             = u2,
    ssrd_MJ_m2_day = 17.1940              # given, per the paper's own note
  )

  out <- Eopen_Penman(df, lat_deg = lat, elev_m = elev, albedo = 0.08,
                       G = 0, wind_fn = "penman1956", verbose = FALSE)

  # Published intermediates (Worked Examples 1-2)
  expect_equal(out$Tmean_C,  11.5,    tolerance = 1e-3)
  expect_equal(out$es,       1.5963,  tolerance = 1e-3)
  expect_equal(out$ea,       0.5614,  tolerance = 1e-3)
  expect_equal(out$vpd,      1.0349,  tolerance = 1e-3)
  expect_equal(out$Ra_MJ,    23.6182, tolerance = 1e-3)
  expect_equal(out$Rso_MJ,   17.9716, tolerance = 1e-3)
  expect_equal(out$Rns_MJ,   15.8184, tolerance = 1e-3)
  expect_equal(out$Rnl_MJ,   7.1784,  tolerance = 5e-3)   # 273.16 vs 273.2, see above
  expect_equal(out$Rn_MJ,    8.6401,  tolerance = 5e-3)

  # Published answer (Worked Example 3). Tolerance set from the observed
  # 0.04% gap, itself attributable to rounding in the paper's own
  # hand-worked intermediate steps.
  expect_equal(out$Eopen_mm, 2.9797, tolerance = 2e-3)
})

test_that("Eopen_Penman default wind function differs from Penman (1956) as documented", {
  # Regression guard on the ~19% gap between the project's default wind
  # function (fao56_shuttleworth) and McMahon's own choice (penman1956)
  # for the same Alice Springs inputs - see the header note in
  # R/Penman_Combination_Open_ETo.R. This is a deliberate methodological
  # choice, not a bug; this test just pins the magnitude so a silent
  # change in either wind function is caught.

  df <- data.frame(
    date           = as.Date("1980-07-20"),
    Tmax_C         = 21.0,
    Tmin_C         = 2.0,
    Tdmean_C       = td_from_ea(0.5614),
    u2             = 0.5903,
    ssrd_MJ_m2_day = 17.1940
  )

  out_default <- Eopen_Penman(df, lat_deg = -23.7951, elev_m = 546,
                               albedo = 0.08, G = 0, verbose = FALSE)
  out_penman  <- Eopen_Penman(df, lat_deg = -23.7951, elev_m = 546,
                               albedo = 0.08, G = 0,
                               wind_fn = "penman1956", verbose = FALSE)

  expect_equal(out_default$Eopen_mm, 3.548, tolerance = 5e-3)
  expect_equal((out_default$Eopen_mm / out_penman$Eopen_mm - 1) * 100,
               19.0, tolerance = 1)
})


# ---------------------------------------------------------------
# 5. Guard behaviour
# ---------------------------------------------------------------

test_that("guards fire on impossible inputs", {
  
  base <- data.frame(
    date           = as.Date("2000-07-01"),
    Tmax_C         = 30, Tmin_C = 18,
    Tdmean_C       = 10, u2 = 2,
    ssrd_MJ_m2_day = 25
  )
  
  bad_rs <- transform(base, ssrd_MJ_m2_day = 500)
  expect_warning(ETo_FAO56(bad_rs, 37.5, 542), "Rs > Ra")
  
  bad_td <- transform(base, Tdmean_C = 35)
  expect_warning(ETo_FAO56(bad_td, 37.5, 542), "Tdew > Tmean")
  
  expect_error(ETo_FAO56(base[, -1], 37.5, 542), "Missing required columns")
  
  expect_error(
    Eopen_Penman(base, 37.5, 542, Tw_C = 20),
    "must be supplied together"
  )
})


# ---------------------------------------------------------------
# 6. ALIGNMENT REGRESSION TEST
# ---------------------------------------------------------------
# Eopen_Penman() sorts df by date internally. G and Tw_C arrive as
# arguments, outside the data frame. If they are not carried through the
# same reordering, the heat storage term is silently attached to the
# wrong days. The failure is invisible: no warning, no NA, plausible
# annual totals, wrong seasonal phasing. Since the seasonal phase lag is
# a headline result, this must stay tested.

make_forcing <- function(n = 60, seed = 1) {
  set.seed(seed)
  dates <- seq(as.Date("2010-01-01"), by = "day", length.out = n)
  J     <- as.integer(strftime(dates, "%j"))
  # Strong seasonality so a misalignment cannot cancel out
  data.frame(
    date           = dates,
    Tmax_C         = 20 + 15 * sin(2 * pi * J / 365),
    Tmin_C         =  8 + 10 * sin(2 * pi * J / 365),
    Tdmean_C       =  2 +  8 * sin(2 * pi * J / 365),
    u2             = 2 + runif(n, -0.5, 0.5),
    ssrd_MJ_m2_day = 18 + 10 * sin(2 * pi * J / 365),
    strd_MJ_m2_day = 25 + 5  * sin(2 * pi * J / 365)
  )
}

test_that("G stays attached to its own date when df arrives unsorted", {
  
  df_sorted <- make_forcing()
  n <- nrow(df_sorted)
  
  # A G series with strong day-to-day structure. A constant or a random
  # series would let a misalignment pass unnoticed.
  G_sorted <- 6 * sin(2 * pi * seq_len(n) / 30)
  
  ref <- Eopen_Penman(df_sorted, 37.5, 542, G = G_sorted, verbose = FALSE)
  
  # Shuffle df and G together, exactly as a user would have them:
  # G[i] belongs to df[i, ], whatever order the rows happen to be in.
  set.seed(99)
  perm <- sample(n)
  shuffled <- Eopen_Penman(df_sorted[perm, ], 37.5, 542,
                           G = G_sorted[perm], verbose = FALSE)
  
  # Both outputs are date-sorted on return, so they should be identical
  expect_equal(shuffled$date,     ref$date)
  expect_equal(shuffled$G_MJ,     ref$G_MJ)
  expect_equal(shuffled$Eopen_mm, ref$Eopen_mm, tolerance = 1e-12)
})

test_that("Tw_C stays attached to its own date when df arrives unsorted", {
  
  df_sorted <- make_forcing()
  n  <- nrow(df_sorted)
  Tw <- 12 + 8 * sin(2 * pi * seq_len(n) / 30)
  
  ref <- Eopen_Penman(df_sorted, 37.5, 542,
                      Tw_C = Tw, Ld_MJ_m2_day = "strd_MJ_m2_day",
                      verbose = FALSE)
  
  set.seed(7)
  perm <- sample(n)
  shuffled <- Eopen_Penman(df_sorted[perm, ], 37.5, 542,
                           Tw_C = Tw[perm], Ld_MJ_m2_day = "strd_MJ_m2_day",
                           verbose = FALSE)
  
  expect_equal(shuffled$Rnl_MJ,   ref$Rnl_MJ,   tolerance = 1e-12)
  expect_equal(shuffled$Eopen_mm, ref$Eopen_mm, tolerance = 1e-12)
})

test_that("the alignment test can actually fail", {
  # Negative control. If a misaligned G produced the same answer as an
  # aligned one, the two tests above would be vacuous. Deliberately
  # break the pairing and confirm the result changes.
  df <- make_forcing()
  n  <- nrow(df)
  G  <- 6 * sin(2 * pi * seq_len(n) / 30)
  
  ok   <- Eopen_Penman(df, 37.5, 542, G = G,       verbose = FALSE)
  bust <- Eopen_Penman(df, 37.5, 542, G = rev(G),  verbose = FALSE)
  
  expect_false(isTRUE(all.equal(ok$Eopen_mm, bust$Eopen_mm)))
})

test_that("wrong-length G and Tw_C are rejected rather than recycled", {
  df <- make_forcing(n = 10)
  expect_error(Eopen_Penman(df, 37.5, 542, G = c(1, 2, 3), verbose = FALSE),
               "length 1 or nrow")
  expect_error(Eopen_Penman(df, 37.5, 542, Tw_C = c(10, 11),
                            Ld_MJ_m2_day = "strd_MJ_m2_day", verbose = FALSE),
               "length 1 or nrow")
})