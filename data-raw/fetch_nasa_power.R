# ===============================================================
# data-raw/fetch_nasa_power.R
#
# Stage C.4 (Ataturk_Master_Checklist_v2.md): NASA POWER as a second,
# independent forcing product for Ataturk Reservoir, 1993-2025.
#
# Checklist is explicit that this is an ENSEMBLE MEMBER, not a
# correction target: POWER (MERRA-2 reanalysis + CERES SYN1deg
# satellite radiation) and ERA5-Land are two independent estimates of
# the same atmospheric forcing. Where they disagree, that disagreement
# IS the result worth reporting (a forcing-uncertainty band), not a bug
# to resolve by picking one. Nothing in this script, or in how its
# output should be used downstream, should make POWER a bias-correction
# input for the ERA5-Land series - they go on the same axis, as two
# lines, per GATE C's own framing.
#
# LOCATION: the reservoir polygon centroid (computed directly from
# `data/ataturk_jrc_max_extent_polygon.geojson` via shapely, 4 Oct
# 2026): lon 38.5878 E, lat 37.6032 N. This is deliberately the
# max-extent-polygon centroid, the SAME geometry already anchoring the
# ERA5-Land re-mask (Stage C.2) and the MODIS LST extraction (Stage
# C.3), so all three products are centred/masked on one consistent
# definition of "the reservoir" rather than each picking its own point.
# (For reference, DAHITI's own target-112 point is lon 38.5648 / lat
# 37.5794 - close but not identical; that's a gauge-point vs.
# polygon-centroid difference, not a conflict to resolve.) POWER's
# native grid is ~0.5 deg, so a single point here already represents
# roughly a 50x50 km cell - a "small grid" of several POWER points
# would mostly resample the same underlying cell and isn't worth the
# extra complexity for an ensemble-member product.
#
# COMMUNITY: "AG" (Agroclimatology) rather than "RE" (Renewable Energy)
# or "SB" (Sustainable Buildings) - AG's parameter set and rounding
# conventions are the ones built for exactly this kind of crop-weather/
# evapotranspiration use case. This is a methods choice, not a
# re-derivable fact; flagged here for Aziz/Tombul to confirm, same as
# the MODIS LST error-threshold choice in fetch_modis_lst.py.
#
# *** THE UNIT PROBLEM, AND HOW THIS SCRIPT HANDLES IT ***
# R/standardise_forcing.R's MAP_NASA_POWER table has carried a warning
# since it was first sketched out: "POWER unit conventions have changed
# across API versions - in particular ALLSKY_SFC_SW_DWN has been served
# in both MJ/m2/day and kW-hr/m2/day, and PS in both kPa and Pa. Do not
# trust this table until you have read the units line of the file you
# downloaded." Checked before writing this script (4 Oct 2026) rather
# than guessing which era's convention applies now:
#   - NASA POWER's own static documentation pages (parameter dictionary,
#     methodology pages) render units through an interactive JS tool,
#     not in static HTML - they could not be read directly.
#   - No outbound network path available to this assistant could reach
#     power.larc.nasa.gov's live API directly to read a real response
#     (blocked by the sandbox's egress policy on both the cloud
#     container and this device) - so the live JSON could not be
#     fetched and inspected ahead of time either.
#   - The one piece of DIRECT evidence available was the `nasapower`
#     package's own vignette, which shows a real `temporal_api="daily"`
#     example output whose decorative header lists
#     CLRSKY_SFC_SW_DWN in **kW-hr/m^2/day** (not MJ/m2/day).
#     ALLSKY_SFC_SW_DWN and ALLSKY_SFC_LW_DWN are the same CERES
#     SYN1deg family of parameters as CLRSKY_SFC_SW_DWN, so the same
#     convention is likely - but "likely" from one sibling parameter's
#     example is not the same as having actually read the units line of
#     THIS download, which is exactly what the original warning asked
#     for. PS was corroborated more directly: pvlib's NASA POWER reader
#     documents explicitly converting "pressure from kPa to Pa" on
#     ingest, i.e. POWER's raw PS is kPa - matching MAP_NASA_POWER's
#     existing factor already.
# Rather than hardcode a guess either way for the radiation unit, this
# script reads the ACTUAL decorative header `get_power()` prints for
# THIS download, extracts the unit string next to each radiation
# parameter with a regex, and converts to the canonical MJ/m2/day only
# based on what that header actually says - loudly, with the raw
# matched string printed, so Aziz sees exactly what was found and
# applied rather than trusting a silent assumption. If the unit string
# can't be confidently classified, the script does NOT guess: it leaves
# the value unconverted and prints an unmissable warning instead, so a
# wrong number is never silently handed downstream into a Penman
# calculation.
#
# RESOLVED 7 Oct 2026: the first real run printed the header, which says
# ALLSKY_SFC_SW_DWN and ALLSKY_SFC_LW_DWN are MJ/m^2/day and PS is kPa, so
# MAP_NASA_POWER in R/standardise_forcing.R was right as written. The
# runtime unit parser described above failed to see the header (it is not
# written to stdout) and was replaced by a magnitude guard.
#
# PREREQUISITE: the `nasapower` package (rOpenSci, CRAN) is not yet in
# this project's renv.lock. Before running this script:
#   install.packages("nasapower")
# and AFTER a successful run, re-snapshot renv with the FULL existing
# package union plus "nasapower" in one explicit call - passing only
# the new package to `packages = c(...)` silently prunes everything
# else not caught by the implicit scanner, exactly as happened with
# `terra` et al. in Stage A (see Ataturk_Master_Checklist_v2.md). The
# known-good explicit list from that episode, plus nasapower:
#   renv::snapshot(packages = c("renv", "here", "rmarkdown", "testthat",
#     "sessioninfo", "dplyr", "tibble", "lubridate", "terra", "tidyr",
#     "purrr", "httr", "jsonlite", "ecmwfr", "keyring", "getPass",
#     "nasapower"))
# Check the printed diff is all `[* -> version]` before confirming, same
# as last time.
# ===============================================================

library(nasapower)

# ---------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------
RESERVOIR_LON <- 38.5878   # JRC max-extent polygon centroid, see header note
RESERVOIR_LAT <- 37.6032

DATE_START <- "1993-01-01"
DATE_END   <- "2025-12-31"

COMMUNITY <- "AG"   # Agroclimatology - see header note; confirm with Tombul

# Canonical mapping needs: T2M_MAX, T2M_MIN, T2MDEW -> temperature/dewpoint
# (standardise_forcing's humidity path, matching ERA5-Land's d2m-based
# approach rather than a direct RH input); PS -> surface pressure;
# ALLSKY_SFC_SW_DWN / ALLSKY_SFC_LW_DWN -> shortwave/longwave radiation;
# PRECTOTCORR -> precipitation; WS10M -> 10 m wind speed, fed through
# standardise_forcing(wind = "from_u10") so it goes through the SAME
# wind10_to_wind2() log-law step as ERA5-Land's u10, rather than using
# POWER's own WS2M (which would use a different, unverified internal
# extrapolation and break the apples-to-apples ensemble comparison).
PARS <- c("T2M_MAX", "T2M_MIN", "T2MDEW", "PS",
          "ALLSKY_SFC_SW_DWN", "ALLSKY_SFC_LW_DWN",
          "PRECTOTCORR", "WS10M")

out_dir  <- here::here("data")
out_path <- file.path(out_dir, "nasa_power_daily_1993-2025.csv")

# ---------------------------------------------------------------
# Fetch (single call - POWER's daily point endpoint is not documented
# with Earth Engine-style request-duration limits; chunk by decade
# manually if this ever errors out on a multi-year span)
# ---------------------------------------------------------------
cat(sprintf("Fetching NASA POWER daily data: community=%s, lon=%.4f, lat=%.4f, %s to %s\n",
            COMMUNITY, RESERVOIR_LON, RESERVOIR_LAT, DATE_START, DATE_END))

power_df <- get_power(
  community    = COMMUNITY,
  lonlat       = c(RESERVOIR_LON, RESERVOIR_LAT),
  pars         = PARS,
  dates        = c(DATE_START, DATE_END),
  temporal_api = "daily"
)

cat(sprintf("Fetched %d rows (%s to %s)\n",
            nrow(power_df), min(power_df$YYYYMMDD), max(power_df$YYYYMMDD)))

# ---------------------------------------------------------------
# Read the units straight off THIS download's own printed header -
# see the long UNIT PROBLEM note above. Do not trust a hardcoded
# assumption for a parameter POWER has changed units for before.
# ---------------------------------------------------------------
# Units VERIFIED against the real header this call prints (first run,
# 7 Oct 2026): T2M_MAX/T2M_MIN/T2MDEW (C), PS (kPa), ALLSKY_SFC_SW_DWN and
# ALLSKY_SFC_LW_DWN (MJ/m^2/day), PRECTOTCORR (mm/day), WS10M (m/s).
# The first version of this script tried to parse those units out of
# capture.output(print(power_df)); that found nothing, because the
# decorative header is not written to stdout, so it could not be captured.
# The CLRSKY kW-hr/m^2/day example in the nasapower vignette did NOT carry
# over to ALLSKY_SFC_SW_DWN here - the real header says MJ/m^2/day, so no
# conversion is needed. Rather than keep a parser that cannot see the
# header, the guard below checks magnitudes: a kW-hr/m^2/day or Pa column
# would sit a factor of 3.6 or 1000 away from these ranges and stop the run.
cat("\n=== Units (verified from the printed header, 7 Oct 2026) ===\n")
cat("  T2M_MAX/T2M_MIN/T2MDEW: C | PS: kPa | ALLSKY_SFC_SW_DWN/LW_DWN: MJ/m^2/day\n")
cat("  PRECTOTCORR: mm/day | WS10M: m/s -> no conversion applied\n")

guard <- function(x, name, lo, hi) {
  m <- mean(x, na.rm = TRUE)
  if (!is.finite(m) || m < lo || m > hi)
    stop(sprintf("%s mean = %.2f is outside the plausible [%g, %g] for this site - ",
                 name, m, lo, hi),
         "POWER's units may have changed; read the printed header before using this file.")
  cat(sprintf("  guard OK: %-18s mean %.2f within [%g, %g]\n", name, m, lo, hi))
}
cat("\n=== Magnitude guard (catches a silent unit change) ===\n")
guard(power_df$ALLSKY_SFC_SW_DWN, "ALLSKY_SFC_SW_DWN", 10, 25)
guard(power_df$ALLSKY_SFC_LW_DWN, "ALLSKY_SFC_LW_DWN", 20, 40)
guard(power_df$PS,                "PS",                80, 105)
guard(power_df$T2M_MAX,           "T2M_MAX",           5, 35)

sw_result <- list(value = power_df$ALLSKY_SFC_SW_DWN)
lw_result <- list(value = power_df$ALLSKY_SFC_LW_DWN)

# ---------------------------------------------------------------
# Build the output frame - raw POWER column names preserved (so this
# file's provenance is traceable), but ALLSKY_SFC_SW_DWN/LW_DWN are now
# GUARANTEED MJ/m2/day (converted above if needed), matching every
# other forcing source's convention in CANONICAL_SCHEMA.
# ---------------------------------------------------------------
out <- data.frame(
  date               = as.Date(power_df$YYYYMMDD),
  T2M_MAX            = power_df$T2M_MAX,
  T2M_MIN            = power_df$T2M_MIN,
  T2MDEW             = power_df$T2MDEW,
  PS                 = power_df$PS,
  ALLSKY_SFC_SW_DWN  = sw_result$value,
  ALLSKY_SFC_LW_DWN  = lw_result$value,
  PRECTOTCORR        = power_df$PRECTOTCORR,
  WS10M              = power_df$WS10M
)

# --- lightweight sanity checks, same spirit as validate_forcing() ---
cat("\n=== Column ranges (eyeball check before trusting this file) ===\n")
for (col in setdiff(names(out), "date")) {
  v <- out[[col]]
  cat(sprintf("  %-20s min=%8.2f  mean=%8.2f  max=%8.2f  NA=%d\n",
              col, min(v, na.rm = TRUE), mean(v, na.rm = TRUE),
              max(v, na.rm = TRUE), sum(is.na(v))))
}

n_dup  <- sum(duplicated(out$date))
gaps   <- as.integer(diff(sort(out$date)))
n_gaps <- sum(gaps > 1)
cat(sprintf("\nDuplicated dates: %d. Gaps >1 day: %d (largest %d day(s)).\n",
            n_dup, n_gaps, if (length(gaps)) max(gaps) else 0L))

write.csv(out, out_path, row.names = FALSE)
cat(sprintf("\nWrote %d rows to %s\n", nrow(out), out_path))
cat("\nReminder: this is an INDEPENDENT ENSEMBLE MEMBER (checklist C.4) - plot it\n")
cat("alongside ERA5-Land, do not use it to bias-correct or replace ERA5-Land values.\n")
