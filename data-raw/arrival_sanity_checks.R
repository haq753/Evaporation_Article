# ===============================================================
# data-raw/arrival_sanity_checks.R
#
# Stage C.1 "Arrival sanity checks" (Ataturk_Master_Checklist_v2.md).
# Runs the 10 checklist bullets against data/daily_forcing_*.csv,
# produced by build_daily_forcing.R. All must pass (or be explained)
# before this data is trusted for Stage D onward.
#
# Reuses the already-verified physics functions (R/utils_physics.R,
# Gate A) rather than reimplementing Ra/Rso/es here - those are tested
# against FAO-56's own worked examples; no reason to duplicate and
# risk drifting from them.
#
# SITE: 37.5 N (Ataturk box centroid-ish), elevation ~540 m - matches
# the bracket already in tests/test_physics.R (41.779 / 15.069
# MJ m^-2 d^-1 at summer/winter solstice).
# ===============================================================

library(dplyr)

source(here::here("R", "utils_physics.R"))
source(here::here("R", "Wind_10m_to_2m.R"))

LAT_DEG <- 37.5
ELEV_M  <- 540

# REVISION (30 Sep 2026): now that build_daily_forcing.R has been run
# for the full 1993:2025 range, picking the file up by name would go
# stale the moment a new one appears. Instead scan data/ for every
# daily_forcing_<start>-<end>.csv and use whichever spans the most
# years - works unmodified whether you just ran the 1993-only test or
# the full 33-year build.
candidates <- list.files(here::here("data"),
                          pattern = "^daily_forcing_[0-9]{4}-[0-9]{4}\\.csv$",
                          full.names = TRUE)
if (length(candidates) == 0) stop("No daily_forcing_*.csv found in data/ - run build_daily_forcing.R first.")
yr <- regmatches(basename(candidates), regexpr("[0-9]{4}-[0-9]{4}", basename(candidates)))
span <- sapply(strsplit(yr, "-"), function(x) as.integer(x[2]) - as.integer(x[1]))
csv_path <- candidates[which.max(span)]

df <- read.csv(csv_path) |> mutate(date = as.Date(date))

cat("Loaded", nrow(df), "days from", basename(csv_path), "\n\n")

# --- expected NA check (30 Sep 2026) -----------------------------
# build_daily_forcing.R's own comment documents this: the accumulated
# variables (ssrd, strd, str, tp, e, pev) need the FOLLOWING day's
# 00:00 value, and the very last calendar day of the whole record has
# no "next" file to pull that from (there is no 2026 forcing file) -
# that day is deliberately left NA rather than dropped or guessed.
# Across a full 1993:2025 run that should be EXACTLY one day: 31 Dec
# 2025. Checking that explicitly here rather than assuming it, since
# "some NAs showed up" could just as easily mean something else broke.
na_rs_dates <- df$date[is.na(df$ssrd_dayTotal)]
cat("Accumulated-variable NA check: ssrd_dayTotal is NA on", length(na_rs_dates),
    "day(s)", if (length(na_rs_dates) > 0) paste0("(", paste(na_rs_dates, collapse = ", "), ")") else "", "\n")
if (length(na_rs_dates) == 1 && na_rs_dates == max(df$date)) {
  cat("  -> matches the expected end-of-record gap (last day, no next-day 00:00 to pull). Fine.\n\n")
} else if (length(na_rs_dates) > 0) {
  cat("  -> DOES NOT match the expected single end-of-record day - investigate before trusting\n")
  cat("     the checks below; na.rm is masking these rather than fixing the cause.\n\n")
} else {
  cat("\n")
}

# --- derived quantities needed for the checks --------------------
df <- df |>
  mutate(
    doy   = as.integer(format(date, "%j")),
    Ra    = Ra_calc(LAT_DEG, doy),
    Rso   = Rso_calc(Ra, ELEV_M),
    Rs    = ssrd_dayTotal,                       # already MJ/m2/day
    Tmax  = t2m_max_val,
    Tmin  = t2m_min_val,
    Tmean = t2m_mean_val,
    Tdmean = d2m_mean_val,
    ea    = es_kPa(Tdmean),                      # ea = es(Tdew), FAO-56 Eq. 14
    es_mean = (es_kPa(Tmax) + es_kPa(Tmin)) / 2,
    VPD   = es_mean - ea,
    # REVISION (29 Sep 2026): run #1 used sqrt(mean(u10)^2+mean(v10)^2)
    # here, i.e. the magnitude of the daily-mean components. That gave
    # a 1.05 m/s mean against the second the checklist's whole 1.5-3
    # expected range - not an isolated bad day, so build_daily_forcing.R
    # was changed to compute sqrt(u10^2+v10^2) at every HOUR first, then
    # daily-mean/max/min THAT (column wspd10_*_val below). This script
    # now just reads that corrected column instead of re-deriving it
    # from components - re-run build_daily_forcing.R first to get a CSV
    # that has wspd10_mean_val before running this.
    u2     = wind10_to_wind2(wspd10_mean_val)
  )

is_summer <- format(df$date, "%m") %in% c("06", "07", "08")
is_winter <- format(df$date, "%m") %in% c("12", "01", "02")

pass <- function(label, ok, detail = "") {
  cat(sprintf("[%s] %s %s\n", if (ok) "PASS" else "CHECK", label, detail))
}

cat("=== Stage C.1 arrival sanity checks (", nrow(df), "days, ", basename(csv_path), ") ===\n\n", sep = "")

# 1. Summer Rs ~27-30, winter Rs ~7-9
cat("1. Solar radiation magnitude\n")
cat("   Summer Rs range:", round(range(df$Rs[is_summer], na.rm = TRUE), 2), "MJ/m2/day (expect ~27-30)\n")
cat("   Winter Rs range:", round(range(df$Rs[is_winter], na.rm = TRUE), 2), "MJ/m2/day (expect ~7-9)\n\n")

# 2. Rs <= Ra every day
# na.rm=TRUE: the one documented end-of-record NA (see check above)
# must not silently turn this whole test result into NA.
n_violate <- sum(df$Rs > df$Ra, na.rm = TRUE)
pass("2. Rs <= Ra on every day", n_violate == 0,
     sprintf("(%d violation(s)%s)", n_violate,
             if (n_violate > 0) paste0(" - worst: ", round(max(df$Rs - df$Ra, na.rm = TRUE), 3), " MJ/m2 over") else ""))

# 3. Rs/Rso in [0.3, 1.0]
ratio <- df$Rs / df$Rso
n_out <- sum(ratio < 0.3 | ratio > 1.0, na.rm = TRUE)
pass("3. Rs/Rso in [0.3, 1.0]", n_out == 0,
     sprintf("(%d day(s) outside range; range found: %.3f - %.3f)",
             n_out, min(ratio, na.rm = TRUE), max(ratio, na.rm = TRUE)))

# 4. Summer Tmax ~38-42, winter Tmin ~0-4
cat("\n4. Temperature extremes\n")
cat("   Summer Tmax range:", round(range(df$Tmax[is_summer], na.rm = TRUE), 2), "C (expect ~38-42)\n")
cat("   Winter Tmin range:", round(range(df$Tmin[is_winter], na.rm = TRUE), 2), "C (expect ~0-4)\n\n")

# 5. Td <= Tmean every day
n_td <- sum(df$Tdmean > df$Tmean, na.rm = TRUE)
pass("5. Td <= Tmean on every day", n_td == 0,
     sprintf("(%d violation(s)%s)", n_td,
             if (n_td > 0) paste0(" - worst: ", round(max(df$Tdmean - df$Tmean, na.rm=TRUE), 2), " C over") else ""))

# 6. Summer VPD ~3-4 kPa
cat("6. Summer VPD range:", round(range(df$VPD[is_summer], na.rm = TRUE), 2), "kPa (expect ~3-4)\n\n")

# 7. u2 mean ~1.5-3 m/s
cat("7. u2 mean:", round(mean(df$u2, na.rm = TRUE), 2), "m/s (expect ~1.5-3); range:",
    round(range(df$u2, na.rm = TRUE), 2), "\n\n")

# 8. Ra bracket at site: summer 41.779, winter 15.069
# REVISION (30 Sep 2026): with multiple years loaded there's one 21 Jun
# and one 21 Dec per year, not a single value - doy for a fixed
# calendar date shifts by 1 in leap years (Ra depends on doy), so the
# per-year values spread by a small, expected amount rather than being
# bit-identical. Check every year's value is within tolerance instead
# of comparing a vector to a scalar (which silently produced the wrong
# answer under all.equal for >1 value).
ra_summer <- df$Ra[format(df$date, "%m-%d") == "06-21"]
ra_winter <- df$Ra[format(df$date, "%m-%d") == "12-21"]
pass("8. Ra bracket (21 Jun = 41.779)", all(abs(ra_summer - 41.779) < 0.05),
     sprintf("(%d year(s), range %.3f-%.3f)", length(ra_summer), min(ra_summer), max(ra_summer)))
pass("   Ra bracket (21 Dec = 15.069)", all(abs(ra_winter - 15.069) < 0.05),
     sprintf("(%d year(s), range %.3f-%.3f)", length(ra_winter), min(ra_winter), max(ra_winter)))

# 9. lake_ice_depth ~0 year-round
# Values of order 1e-18 or smaller are floating-point noise around a
# true zero, not ice - only count magnitudes above that as real ice.
licd_nonzero <- df$licd_mean_val > 1e-6 & !is.na(df$licd_mean_val)
cat("9. licd_mean_val range:", round(range(df$licd_mean_val, na.rm = TRUE), 4), "m (expect ~0 year-round)\n")
cat("   Non-noise (>1e-6 m) on", sum(licd_nonzero), "of", nrow(df), "days",
    if (sum(licd_nonzero) > 0)
      sprintf(" (%s to %s)", min(df$date[licd_nonzero]), max(df$date[licd_nonzero]))
    else "", "\n\n")

# 10. lake_mix_layer_depth - check for abrupt jumps
lmld_diff <- diff(df$lmld_mean_val)
big_jump <- which(abs(lmld_diff) > 2)   # >2 m day-to-day jump = suspicious
cat("10. lmld_mean_val day-to-day jumps >2m:", length(big_jump), "occurrence(s)\n")
if (length(big_jump) > 0) {
  show <- big_jump[seq_len(min(30, length(big_jump)))]
  cat("    at:", paste(as.character(df$date[show + 1]), collapse = ", "),
      if (length(big_jump) > 30) sprintf(" ... +%d more", length(big_jump) - 30) else "", "\n")
}

cat("\n=== Done. PASS items are solid; CHECK items need a look before Stage D. ===\n")
cat("Numeric-range items (1, 4, 6, 7, 9) have no pass/fail threshold coded -\n")
cat("read them against the expected ranges above and use judgement, especially\n")
cat("given the box-mean-vs-polygon caveat already flagged in build_daily_forcing.R.\n")

# ===============================================================
# EXTRA DIAGNOSTIC (rev 2, 29 Sep 2026): run #1's version of this
# block used "Tmin < 0 C" as its cold threshold, which flagged ~1/3
# of the year (ordinary winter nights below freezing are NOT
# anomalous for this site - Adiyaman, inside the request box, has a
# documented Dec-Feb mean daily minimum of +1.8 to +3.6 C, so plain
# sub-zero nights are unremarkable and drowned out the days that
# actually matter). Tightened to Tmin < -7 C - clearly colder than
# that normal range - to isolate genuine cold-snap days instead of
# printing a third of the year.
# ===============================================================
cat("\n=== Extra diagnostic: do the flagged days line up? ===\n\n")

# REVISION (1 Oct 2026): the full 33-year run surfaced a new question
# the 1993-only test couldn't: 465 days fall outside Rs/Rso [0.3, 1.0],
# some down to 0.066 - well below what FAO-56 treats as "heavily
# overcast" (~0.3). The 15 lowest-ratio days are mostly NOT cold or
# icy (several sit at Tmin 5-10 C with licd at FP-noise levels), so
# they are not the same winter cold-snap mechanism already explained.
# Two live hypotheses, not yet distinguished: (a) genuine heavy-rain/
# storm days (thick rain cloud can plausibly cut Rs this much), or
# (b) something else - e.g. dust/haze events, which are a real,
# documented phenomenon in this region (Mesopotamian dust storms) but
# would need a different data source to confirm since ERA5-Land's
# radiation here is model-derived, not a dust measurement. Added
# tp_dayTotal (precipitation, mm) to this table so the two can be told
# apart from data already on hand: a low-ratio day with high tp is a
# real storm (expected); a low-ratio day with near-zero tp needs more
# digging before calling it explained.
flagged <- df |>
  mutate(ratio = Rs / Rso) |>
  filter(ratio < 0.3 | Tmin < -7 | licd_mean_val > 1e-6) |>
  select(date, Rs, Ra, Rso, ratio, Tmin, Tmax, VPD, licd_mean_val, tp_dayTotal) |>
  arrange(date)

cat("Days with Rs/Rso < 0.3, or Tmin < -7 C, or non-noise licd (",
    nrow(flagged), "day(s)):\n\n", sep = "")

# REVISION (30 Sep 2026): the full multi-year file can flag hundreds of
# days - dumping every row was only readable for the 1993 single-year
# test. Full row-by-row printout stays for a short run; a long run
# gets a per-year count plus the most extreme rows instead, which is
# what you actually need to judge "is this normal winter spread or a
# real problem" without scrolling hundreds of lines.
if (nrow(df) <= 400) {
  print(as.data.frame(flagged), row.names = FALSE)
} else {
  cat("(more than one year loaded - showing a summary, not every row)\n\n")
  by_year <- flagged |>
    mutate(year = format(date, "%Y")) |>
    count(year, name = "n_flagged_days")
  cat("Flagged days per year:\n")
  print(as.data.frame(by_year), row.names = FALSE)

  cat("\n15 lowest Rs/Rso days (most overcast relative to clear sky):\n")
  print(as.data.frame(flagged |> arrange(ratio) |> head(15)), row.names = FALSE)

  cat("\n15 coldest Tmin days:\n")
  print(as.data.frame(flagged |> arrange(Tmin) |> head(15)), row.names = FALSE)
}

cat("\nRead this as: an isolated overcast day (low ratio, ordinary Tmin,\n")
cat("no ice) is just weather. A cluster where LOW ratio, very LOW Tmin,\n")
cat("AND non-zero licd land on the SAME or adjacent dates is one\n")
cat("coherent cold-snap episode - expected, not a bug. Rows that show\n")
cat("one extreme in isolation, with no support from the others nearby,\n")
cat("are the ones worth chasing individually. A year with a much higher\n")
cat("n_flagged_days than its neighbours (in the per-year table above) is\n")
cat("worth a closer look on its own.\n")

# ===============================================================
# SUMMER Tmax/VPD LOW TAIL (1 Oct 2026): left unresolved in the 1993-
# only run - summer Tmax as low as 21.29 C and summer VPD as low as
# 0.65 kPa, neither checked against anything at the time. Same move as
# the winter Rs/Rso tail above: look at tp_dayTotal on the actual
# lowest days instead of assuming "probably a rainy summer day" is the
# answer without checking.
# ===============================================================
cat("\n=== Summer Tmax/VPD low tail (flagged unresolved in the 1993-only test) ===\n\n")

summer_extremes <- df |>
  filter(is_summer) |>
  mutate(ratio = Rs / Rso) |>
  select(date, Tmax, Tmin, Rs, ratio, VPD, tp_dayTotal)

cat("5 lowest summer Tmax days:\n")
print(as.data.frame(summer_extremes |> arrange(Tmax) |> head(5)), row.names = FALSE)

cat("\n5 lowest summer VPD days:\n")
print(as.data.frame(summer_extremes |> arrange(VPD) |> head(5)), row.names = FALSE)

cat("\nSame read as before: low Tmax/VPD together with high tp_dayTotal and\n")
cat("a low ratio is a rainy/overcast summer day (expected, same mechanism\n")
cat("as the winter storms above). Low Tmax/VPD with near-zero rain and a\n")
cat("normal ratio would NOT fit that story and needs a closer look.\n")
