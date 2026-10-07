# ===============================================================
# data-raw/build_daily_forcing.R
#
# Turns the 792 GRIB files in data-raw/era5land/ into one daily data
# frame: date, Tmax_C, Tmin_C, Tmean_C, Tdmean_C, u2 (2m wind, m/s),
# ssrd_MJ_m2_day, strd_MJ_m2_day, str_MJ_m2_day, tp_mm_day, e_mm_day,
# pev_mm_day, sp_kPa, skt_C, plus the 7 lake state variables daily-
# meaned the same way. This is the input `prepare_forcing()` in
# R/utils_physics.R needs (date, Tmax_C, Tmin_C, Tdmean_C, u2,
# ssrd_MJ_m2_day) - check that file's exact expected column names
# before wiring this into ETo_FAO56()/Eopen_Penman(), they aren't
# duplicated here.
#
# Everything below is built from THREE things confirmed against real
# output this session, not assumed - each one would have silently
# produced wrong numbers if guessed instead of checked:
#
# ---------------------------------------------------------------
# 1. LAYER ORDER - confirmed via data-raw/inspect_grib_structure.R,
#    29 Sep 2026, on era5land_forcing_1993_01.grib and
#    era5land_lake_1993_01.grib: layers are INTERLEAVED BY TIME (layer
#    1 = variable A hour 0, layer 2 = variable B hour 0, ... layer 12 =
#    variable L hour 0, layer 13 = variable A hour 1, ...), not grouped
#    by variable. Every variable's layer indices came back non-
#    contiguous, confirming this. So variables are identified here by
#    NAME (varnames(), not position, and matched against `time()` for
#    the timestamp - never by "layer N belongs to variable N/12".
#
# 2. VARIABLE NAMES - GDAL's GRIB driver gives clean shortName-based
#    labels (e.g. "2T_0-SFC") for parameters in its built-in tables,
#    but ECMWF's own local-table-228 parameters (all 7 lake variables,
#    plus potential_evaporation) fall through to a raw
#    "varN of table 228 of center ECMWF" label instead - confirmed from
#    the same inspection run. The N in that label IS the last 1-3
#    digits of the ECMWF paramId (e.g. "var251" = paramId 228251).
#    Cross-checked against ECMWF's own parameter database
#    (confluence.ecmwf.int, ERA5-Land data documentation, Table 2/3),
#    not inferred from request order - request order would have been
#    WRONG here (lake_bottom_temperature and lake_total_layer_temperature
#    are adjacent paramIds but in the opposite order from how they were
#    requested; lake_shape_factor sits between them and lake_ice_*).
#
# 3. ACCUMULATION CONVENTION - confirmed against ECMWF's ERA5-Land
#    documentation (not plain ERA5, which works differently): ssrd,
#    str, strd, tp, e and pev are running totals since 00 UTC of that
#    day, not per-hour deltas. The value timestamped 00:00 on day D+1
#    IS day D's full total - summing 24 hourly values instead would
#    massively over-count (this is exactly what the master checklist's
#    own sanity check anticipated: "summer Rs ~300 means running
#    accumulations were summed"). The 7 lake variables, by contrast,
#    are confirmed INSTANTANEOUS (ECMWF's Table 2) - daily-mean them
#    normally, no end-of-day lookup needed.
#    Source: https://confluence.ecmwf.int/display/CKB/ERA5-Land%3A+data+documentation
#
# ---------------------------------------------------------------
# STATUS: written from verified structure, NOT yet run against the
# real 33-year archive (this session has no R/terra environment to
# test in). Run it on 1993 ALONE first (see RUN_YEARS below) and sanity
# -check that one year's output before trusting the full 1993-2025
# loop - cheap to check now, expensive to discover wrong after a full
# run.
#
# REVISION 1 (29 Sep 2026): 1993 test run completed, arrival_sanity_
# checks.R run against it. Fixed the u2 (wind speed) aggregation - see
# the note inside build_month() below - after the sanity check's u2
# mean (1.05 m/s) came back below the entire 1.5-3 m/s expected range,
# for the whole year, not just isolated days. Re-run 1993 again after
# this change and re-check item 7 before widening RUN_YEARS.
# ===============================================================

library(terra)
library(dplyr)
library(tidyr)
library(purrr)
library(lubridate)

era_dir <- here::here("data-raw", "era5land")
out_dir <- here::here("data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- verified variable lookup table -----------------------------
# base_label = varnames(r) with terra's trailing "_<n>" duplicate-layer
# suffix stripped (sub("_[0-9]+$", "", varnames(r))).
VAR_TABLE <- tribble(
  ~base_label,                                  ~code,   ~kind,           ~unit_in, ~unit_out,
  "2T_0-SFC",                                    "t2m",   "instantaneous", "K",      "C",
  "2D_0-SFC",                                    "d2m",   "instantaneous", "K",      "C",
  "SKT_0-SFC",                                   "skt",   "instantaneous", "K",      "C",
  "SP_0-SFC",                                    "sp",    "instantaneous", "Pa",     "kPa",
  "10U_0-SFC",                                   "u10",   "instantaneous", "m/s",    "m/s",
  "10V_0-SFC",                                   "v10",   "instantaneous", "m/s",    "m/s",
  "SSRD_0-SFC",                                  "ssrd",  "accumulated",   "J/m2",   "MJ/m2",
  "STR_0-SFC",                                   "str",   "accumulated",   "J/m2",   "MJ/m2",
  "STRD_0-SFC",                                  "strd",  "accumulated",   "J/m2",   "MJ/m2",
  "TP_0-SFC",                                    "tp",    "accumulated",   "m",      "mm",
  "E_0-SFC",                                     "e",     "accumulated",   "m",      "mm",
  "var251 of table 228 of center ECMWF_0-SFC",   "pev",   "accumulated",   "m",      "mm",
  "var8 of table 228 of center ECMWF_0-SFC",     "lmlt",  "instantaneous", "K",      "C",
  "var9 of table 228 of center ECMWF_0-SFC",     "lmld",  "instantaneous", "m",      "m",
  "var10 of table 228 of center ECMWF_0-SFC",    "lblt",  "instantaneous", "K",      "C",
  "var11 of table 228 of center ECMWF_0-SFC",    "ltlt",  "instantaneous", "K",      "C",
  "var12 of table 228 of center ECMWF_0-SFC",    "lshf",  "instantaneous", "-",      "-",
  "var13 of table 228 of center ECMWF_0-SFC",    "lict",  "instantaneous", "K",      "C",
  "var14 of table 228 of center ECMWF_0-SFC",    "licd",  "instantaneous", "m",      "m",
)

# NOTE ON THE "_0-SFC" SUFFIX: seen consistently in the inspection run
# for every variable. If a future terra/GDAL version on your machine
# prints these labels differently, `extract_month()` below will throw
# "unrecognized variable label" rather than silently mis-assigning -
# update VAR_TABLE's base_label column to match, don't suppress the error.

strip_suffix <- function(x) sub("_[0-9]+$", "", x)

# --- load one (year, month, tag) GRIB file into a long data frame ----
# Returns: label, code, kind, time (POSIXct), value (converted to
# unit_out), one row per (variable, hour).
load_grib_long <- function(year, month, tag) {
  path <- file.path(era_dir, sprintf("era5land_%s_%d_%02d.grib", tag, year, month))
  if (!file.exists(path)) stop("Missing file: ", path)

  r <- rast(path)
  base_lab <- strip_suffix(varnames(r))
  tm <- time(r)

  unmatched <- setdiff(unique(base_lab), VAR_TABLE$base_label)
  if (length(unmatched) > 0) {
    stop("unrecognized variable label(s) in ", basename(path), ": ",
         paste(unmatched, collapse = ", "),
         " - update VAR_TABLE, do not guess the mapping.")
  }

  # spatial mean per layer (box mean - Stage C.2 will remask to the
  # actual reservoir polygon; this is a placeholder consistent with
  # "reservoir-mean forcing, not a point" per the checklist, pending
  # the real polygon)
  vals <- global(r, "mean", na.rm = TRUE)[[1]]

  df <- tibble(base_label = base_lab, time = tm, value_raw = vals) |>
    left_join(VAR_TABLE, by = "base_label")

  # unit conversion
  df <- df |>
    mutate(value = case_when(
      unit_in == "K" & unit_out == "C"   ~ value_raw - 273.15,
      unit_in == "Pa" & unit_out == "kPa" ~ value_raw / 1000,
      unit_in == "J/m2" & unit_out == "MJ/m2" ~ value_raw / 1e6,
      unit_in == "m" & unit_out == "mm"  ~ value_raw * 1000,
      TRUE ~ value_raw
    ))

  # sanity: expected layer count per variable = hours in this month
  n_hours_expected <- length(unique(tm))
  bad <- df |> count(code) |> filter(n != n_hours_expected)
  if (nrow(bad) > 0) {
    stop("layer-count mismatch in ", basename(path), " for: ",
         paste(bad$code, collapse = ", "),
         " - expected ", n_hours_expected, " hourly layers each.")
  }

  df |> select(code, kind, time, value)
}

# --- accumulated-variable end-of-day extraction ------------------
# For each calendar day in `long_df`, the correct daily total is the
# value timestamped 00:00 on the NEXT day (see REVISION note above,
# point 3). `next_day_00h` is a named list/vector of {date: value} for
# the accumulated codes, built by peeking at the following month's
# first timestep (or the following year's Jan file, at a Dec/Jan
# boundary). The LAST calendar day of the whole record (31 Dec 2025)
# has no "1 Jan 2026" file to pull from - handled by leaving it NA and
# printing a warning, not by silently dropping or estimating it.
get_next_day_00h <- function(year, month) {
  ny <- year; nm <- month + 1
  if (nm > 12) { ny <- year + 1; nm <- 1 }
  path <- file.path(era_dir, sprintf("era5land_forcing_%d_%02d.grib", ny, nm))
  if (!file.exists(path)) return(NULL)   # end of record - handled by caller
  df <- load_grib_long(ny, nm, "forcing") |> filter(kind == "accumulated")
  df |> filter(time == min(time)) |> select(code, value) |> rename(next00 = value)
}

# --- build one month's daily rows --------------------------------
build_month <- function(year, month) {
  forcing <- load_grib_long(year, month, "forcing")
  lake    <- load_grib_long(year, month, "lake")

  # REVISION (29 Sep 2026, after arrival_sanity_checks.R run #1):
  # u2 must be built from PER-HOUR wind SPEED, then daily-aggregated -
  # NOT from the magnitude of the already daily-meaned u10/v10
  # components. Vector-averaging components first underestimates true
  # mean speed whenever wind direction varies within the day (it
  # allows opposing-direction hours to partially cancel before the
  # magnitude is even taken). This was flagged as a risk when
  # arrival_sanity_checks.R was written, and check #7 confirmed it:
  # components method gave a 1.05 m/s mean against a 1.5-3 m/s
  # expected floor - not one bad day, the whole-year mean, which is
  # the signature of a systematic method bias rather than real calm
  # weather. Fix: compute sqrt(u10^2+v10^2) at EVERY hour first, then
  # let that hourly series go through the same daily mean/max/min
  # aggregation as every other instantaneous variable.
  wind_hourly <- forcing |>
    filter(code %in% c("u10", "v10")) |>
    select(code, time, value) |>
    pivot_wider(names_from = code, values_from = value) |>
    mutate(code = "wspd10", kind = "instantaneous",
           value = sqrt(u10^2 + v10^2)) |>
    select(code, kind, time, value)

  inst <- bind_rows(forcing, lake, wind_hourly) |>
    filter(kind == "instantaneous") |>
    mutate(date = as.Date(time)) |>
    group_by(date, code) |>
    summarise(
      mean_val = mean(value, na.rm = TRUE),
      max_val  = max(value, na.rm = TRUE),
      min_val  = min(value, na.rm = TRUE),
      .groups = "drop"
    )

  # daily instantaneous summary, wide
  daily_inst <- inst |>
    pivot_wider(id_cols = date,
                names_from = code,
                values_from = c(mean_val, max_val, min_val),
                names_glue = "{code}_{.value}")

  # accumulated variables: pull the day-D+1 00:00 value for each day D
  next00 <- get_next_day_00h(year, month)
  acc <- forcing |> filter(kind == "accumulated")
  days_in_month <- unique(as.Date(acc$time))

  acc_daily <- purrr::map_dfr(days_in_month, function(d) {
    next_d <- d + 1
    if (next_d %in% as.Date(acc$time)) {
      row <- acc |> filter(as.Date(time) == next_d, format(time, "%H:%M") == "00:00")
    } else if (!is.null(next00)) {
      row <- next00 |> mutate(date = d) |> select(code, value = next00, date)
      return(row |> pivot_wider(id_cols = date, names_from = code,
                                 values_from = value, names_glue = "{code}_dayTotal"))
    } else {
      return(tibble(date = d))   # end of record, no next-day value available
    }
    row |> mutate(date = d) |>
      select(date, code, value) |>
      pivot_wider(id_cols = date, names_from = code, values_from = value,
                  names_glue = "{code}_dayTotal")
  })

  full_join(daily_inst, acc_daily, by = "date")
}

# ---------------------------------------------------------------
# RUN_YEARS - test on 1993 alone first. Widen to 1993:2025 only after
# checking January/July 1993 by hand against the checklist's arrival
# sanity-check ranges (Stage C.1, Ataturk_Master_Checklist_v2.md).
# ---------------------------------------------------------------
RUN_YEARS <- 1993:2025

all_months <- expand.grid(month = 1:12, year = RUN_YEARS)
result <- purrr::map2_dfr(all_months$year, all_months$month, build_month) |>
  arrange(date)

out_path <- file.path(out_dir, sprintf("daily_forcing_%s.csv",
                                        paste(range(RUN_YEARS), collapse = "-")))
write.csv(result, out_path, row.names = FALSE)
cat("Wrote", nrow(result), "daily rows to", out_path, "\n")
cat("\nColumn names:\n"); print(names(result))
cat("\nFirst few rows:\n"); print(head(result))
cat("\nQuick range check (compare against Stage C.1 sanity ranges):\n")
cat("t2m_max_val range:", range(result$t2m_max_val, na.rm = TRUE), "C\n")
cat("t2m_min_val range:", range(result$t2m_min_val, na.rm = TRUE), "C\n")
cat("ssrd_dayTotal range:", range(result$ssrd_dayTotal, na.rm = TRUE), "MJ/m2\n")
