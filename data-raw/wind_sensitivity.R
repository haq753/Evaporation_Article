# ===============================================================
# data-raw/wind_sensitivity.R   (8 Oct 2026)
#
# Question: how much does the ERA5-Land vs NASA POWER 10 m wind gap
# (POWER ~1.77x higher, r = 0.84) move open-water evaporation?
#
# Everything except wind is held fixed at the ERA5-Land reservoir-mean
# forcing, so the E difference is attributable to wind alone. Two wind
# functions are crossed with the wind source, because the fao56_shuttleworth
# vs penman1956 choice is itself a flagged open decision.
#
# Reads only local files. Writes data/wind_sensitivity_results.csv.
# Run from the RStudio project:  source("data-raw/wind_sensitivity.R")
# ===============================================================
library(dplyr)

source(here::here("R", "utils_physics.R"))
source(here::here("R", "Wind_10m_to_2m.R"))
source(here::here("R", "standardise_forcing.R"))
source(here::here("R", "Penman_Combination_Open_ETo.R"))

LAT_DEG <- 37.6032   # JRC polygon centroid (Main.qmd uses 37.5; see checklist)
ELEV_M  <- 542       # as Main.qmd
GLEV_MM_YR_1993_2018 <- 1292   # GLEV, Hylak_id 1348, from C.5

# ---- ERA5-Land reservoir-mean, canonical schema ------------------------
# The reservoir file is already in deg C / kPa / MJ m-2 day-1 (verified in
# C.4), so no unit conversion here except kPa -> Pa for sp.
era <- read.csv(here::here("data", "daily_forcing_reservoir_1993-2025.csv"))
era$date <- as.Date(era$date)
era_can <- data.frame(
  date           = era$date,
  Tmax_C         = era$t2m_max_val,
  Tmin_C         = era$t2m_min_val,
  Tdmean_C       = era$d2m_mean_val,
  sp_Pa          = era$sp_mean_val * 1000,
  ssrd_MJ_m2_day = era$ssrd_dayTotal,
  wspd10_era     = era$wspd10_mean_val
)
stopifnot(all(era_can$Tmax_C > -40, na.rm = TRUE), all(era_can$Tmax_C < 60, na.rm = TRUE),
          mean(era_can$sp_Pa, na.rm = TRUE) > 90000,
          mean(era_can$ssrd_MJ_m2_day, na.rm = TRUE) > 10)

pw <- read.csv(here::here("data", "nasa_power_daily_1993-2025.csv"))
pw$date <- as.Date(pw$date)
stopifnot(identical(era_can$date, pw$date))
era_can$wspd10_pow <- pw$WS10M

# ---- one run = (10 m wind column) x (wind function) --------------------
run_E <- function(base, u10, wind_fn, scale = 1) {
  d <- base
  d$u2 <- wind10_to_wind2(u10 * scale, method = "fao56")
  o <- suppressWarnings(suppressMessages(
         Eopen_Penman(d, lat_deg = LAT_DEG, elev_m = ELEV_M,
                      wind_fn = wind_fn, verbose = FALSE)))
  o$year <- as.integer(format(o$date, "%Y"))
  # 2025-12-31 has NA for every accumulated variable (ssrd, strd, ...): the
  # last ERA5-Land accumulation window needs the next day's 00 UTC step, which
  # the download does not include. Annual sums therefore use COMPLETE years only
  # (1993-2024); 2025 is dropped here, not silently zero-filled.
  ann <- o |> group_by(year) |>
    summarise(E_mm = sum(Eopen_mm), n_valid = sum(!is.na(Eopen_mm)), .groups = "drop") |>
    filter(n_valid >= 365)
  stopifnot(nrow(ann) == 32, max(ann$year) == 2024)
  list(ann = ann, u2_mean = mean(d$u2, na.rm = TRUE),
       aero_share = sum(o$aero_mm[o$year <= 2024], na.rm = TRUE) /
                    sum(o$Eopen_mm[o$year <= 2024], na.rm = TRUE))
}

summarise_run <- function(label, r) {
  a  <- r$ann
  data.frame(
    run = label,
    u2_mean_ms = round(r$u2_mean, 3),
    E_mm_yr_1993_2024 = round(mean(a$E_mm), 0),
    E_mm_yr_1993_2018 = round(mean(a$E_mm[a$year <= 2018]), 0),
    aero_share = round(r$aero_share, 3)
  )
}

res <- list()
for (wf in c("fao56_shuttleworth", "penman1956")) {
  res[[paste0("ERA5L wind | ", wf)]]  <- run_E(era_can, era_can$wspd10_era, wf)
  res[[paste0("POWER wind | ", wf)]]  <- run_E(era_can, era_can$wspd10_pow, wf)
}
# Scaling ladder on ERA5-Land wind (fao56_shuttleworth): dE/du
for (s in c(0.75, 1.25, 1.5, 1.75)) {
  res[[sprintf("ERA5L wind x%.2f | fao56_shuttleworth", s)]] <-
    run_E(era_can, era_can$wspd10_era, "fao56_shuttleworth", scale = s)
}

out <- bind_rows(lapply(names(res), function(n) summarise_run(n, res[[n]])))
out$vs_GLEV_1993_2018_mm <- out$E_mm_yr_1993_2018 - GLEV_MM_YR_1993_2018

cat("\nLAT_DEG =", LAT_DEG, " ELEV_M =", ELEV_M,
    " G = 0, FAO-56 Eq.39 longwave branch, albedo", ALBEDO_WATER, "\n")
cat("GLEV 1993-2018 benchmark:", GLEV_MM_YR_1993_2018, "mm/yr\n\n")
print(out, row.names = FALSE)

dir.create(here::here("data"), showWarnings = FALSE)
write.csv(out, here::here("data", "wind_sensitivity_results.csv"), row.names = FALSE)
cat("\nWrote data/wind_sensitivity_results.csv\n")
