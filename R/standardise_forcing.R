# ===============================================================
# R/standardise_forcing.R
#
# One canonical schema. Every data source gets an adapter that maps onto
# it. Applied once at load; nothing downstream ever touches a source
# column name again.
#
# The physics functions are NOT name-agnostic — prepare_forcing() in
# utils_physics.R requires the canonical names below. That is deliberate:
# one schema enforced in one place beats every function guessing.
#
# Adding a source means adding a mapping table here and nothing else.
# ===============================================================

library(dplyr)


# ---------------------------------------------------------------
# 1. Canonical schema
# ---------------------------------------------------------------

CANONICAL_SCHEMA <- tibble::tribble(
  ~name,                ~unit,          ~required, ~plausible_lo, ~plausible_hi,
  "date",               "Date",          TRUE,      NA,            NA,
  "Tmax_C",             "deg C",         TRUE,     -40,            60,
  "Tmin_C",             "deg C",         TRUE,     -40,            50,
  "Tdmean_C",           "deg C",         TRUE,     -50,            40,
  "u2",                 "m/s",           TRUE,       0,            30,
  "ssrd_MJ_m2_day",     "MJ/m2/day",     TRUE,       0,            45,
  "sp_Pa",              "Pa",           FALSE,   50000,        110000,
  "strd_MJ_m2_day",     "MJ/m2/day",    FALSE,       5,            50,
  "Rnl_era5_MJ_m2_day", "MJ/m2/day",    FALSE,      -5,            20,
  "tp_mm",              "mm",           FALSE,       0,           300,
  "e_era5_mm",          "mm",           FALSE,      -5,            25,
  "pev_era5_mm",        "mm",           FALSE,      -5,            30,
  "skt_C",              "deg C",        FALSE,     -40,            70,
  "Tw_C",               "deg C",        FALSE,      -2,            40,
  "d_mix_m",            "m",            FALSE,       0,            60,
  "Ttot_C",             "deg C",        FALSE,      -2,            40,
  "Tbot_C",             "deg C",        FALSE,      -2,            40,
  "ice_m",              "m",            FALSE,       0,             3,
  "lake_cover",         "fraction",     FALSE,       0,             1
)


# ---------------------------------------------------------------
# 2. Source adapters
# ---------------------------------------------------------------
# `factor` and `offset` give:  canonical = source * factor + offset
#
# Sign conventions worth reading twice. ERA5 accumulates downward-positive,
# so `str` (net thermal) and `e` (evaporation) both arrive NEGATIVE. The
# factors below flip them so that outgoing longwave and evaporation are
# positive, matching the convention used in Rnl_water() and everywhere in
# the analysis. Getting this wrong inverts the longwave cross-check.

MAP_ERA5_LAND <- tibble::tribble(
  ~source,          ~name,                 ~factor,  ~offset,
  "date",           "date",                     NA,       NA,   # passthrough
  "t2m_max",        "Tmax_C",                    1,  -273.15,
  "t2m_min",        "Tmin_C",                    1,  -273.15,
  "d2m_mean",       "Tdmean_C",                  1,  -273.15,
  "sp_mean",        "sp_Pa",                     1,        0,
  "ssrd",           "ssrd_MJ_m2_day",         1e-6,        0,
  "strd",           "strd_MJ_m2_day",         1e-6,        0,
  "str",            "Rnl_era5_MJ_m2_day",    -1e-6,        0,   # sign flip
  "tp",             "tp_mm",                  1000,        0,
  "e",              "e_era5_mm",             -1000,        0,   # sign flip
  "pev",            "pev_era5_mm",           -1000,        0,   # sign flip
  "skt_mean",       "skt_C",                     1,  -273.15,
  "lmlt_mean",      "Tw_C",                      1,  -273.15,
  "lmld_mean",      "d_mix_m",                   1,        0,
  "ltlt_mean",      "Ttot_C",                    1,  -273.15,
  "lblt_mean",      "Tbot_C",                    1,  -273.15,
  "licd_mean",      "ice_m",                     1,        0,
  "cl",             "lake_cover",                1,        0
)

# Units VERIFIED 7 Oct 2026 against the header printed by the real POWER
# download (community AG, daily, nasapower::get_power): T2M_MAX / T2M_MIN /
# T2MDEW in C, PS in kPa, ALLSKY_SFC_SW_DWN and ALLSKY_SFC_LW_DWN in
# MJ/m^2/day, PRECTOTCORR in mm/day, WS10M in m/s. POWER has changed unit
# conventions across API versions before, so data-raw/fetch_nasa_power.R
# also carries a magnitude guard that stops the run if they change again.
MAP_NASA_POWER <- tibble::tribble(
  ~source,             ~name,              ~factor, ~offset,
  "date",              "date",                  NA,      NA,
  "T2M_MAX",           "Tmax_C",                 1,       0,
  "T2M_MIN",           "Tmin_C",                 1,       0,
  "T2MDEW",            "Tdmean_C",               1,       0,
  "PS",                "sp_Pa",               1000,       0,   # kPa -> Pa
  "ALLSKY_SFC_SW_DWN", "ssrd_MJ_m2_day",         1,       0,   # MJ/m2/day, verified
  "ALLSKY_SFC_LW_DWN", "strd_MJ_m2_day",         1,       0,
  "PRECTOTCORR",       "tp_mm",                  1,       0
)

SOURCE_MAPS <- list(
  era5_land  = MAP_ERA5_LAND,
  nasa_power = MAP_NASA_POWER
)


# ---------------------------------------------------------------
# 3. The standardiser
# ---------------------------------------------------------------

#' Map a source data frame onto the canonical schema
#'
#' @param df          daily data from one source
#' @param source      "era5_land" or "nasa_power"
#' @param wind        how to obtain u2:
#'                    "from_uv"  — expects u10_mean and v10_mean components
#'                    "from_u10" — expects a single u10 column (name in `u10_col`)
#'                    "present"  — u2 already exists, leave alone
#' @param u10_col     column holding 10 m wind speed when wind = "from_u10"
#' @param wind_method passed to wind10_to_wind2(): "fao56" or "log"
#' @param keep_extra  retain unmapped source columns (default TRUE)
#' @param strict      error rather than warn on implausible ranges
#' @return data frame using canonical names and units
standardise_forcing <- function(df,
                                source = c("era5_land", "nasa_power"),
                                wind = c("from_uv", "from_u10", "present"),
                                u10_col = "u10",
                                wind_method = "fao56",
                                keep_extra = TRUE,
                                strict = FALSE) {
  
  source <- match.arg(source)
  wind   <- match.arg(wind)
  map    <- SOURCE_MAPS[[source]]
  
  if (!exists("wind10_to_wind2") && wind != "present")
    stop("Source R/Wind_10m_to_2m.R before standardising, or use wind = 'present'.")
  
  # --- date -----------------------------------------------------
  if (!"date" %in% names(df))
    stop("`df` needs a `date` column. Rename it before standardising: ",
         "the date is the one thing this function will not guess at.")
  out <- data.frame(date = as.Date(df$date))
  
  # --- mapped variables ----------------------------------------
  applied <- character()
  for (i in seq_len(nrow(map))) {
    src <- map$source[i]
    if (src == "date" || !src %in% names(df)) next
    out[[map$name[i]]] <- df[[src]] * map$factor[i] + map$offset[i]
    applied <- c(applied, src)
  }
  
  # --- wind ----------------------------------------------------
  if (wind == "from_uv") {
    uv <- c("u10_mean", "v10_mean")
    if (!all(uv %in% names(df)))
      stop("wind = 'from_uv' needs columns ", paste(uv, collapse = " and "),
           ". Found: ", paste(intersect(uv, names(df)), collapse = ", "))
    u10     <- sqrt(df$u10_mean^2 + df$v10_mean^2)
    out$u2  <- wind10_to_wind2(u10, method = wind_method)
    out$u10 <- u10
    applied <- c(applied, uv)
  } else if (wind == "from_u10") {
    if (!u10_col %in% names(df))
      stop("Column '", u10_col, "' not found for wind = 'from_u10'.")
    out$u2  <- wind10_to_wind2(df[[u10_col]], method = wind_method)
    out$u10 <- df[[u10_col]]
    applied <- c(applied, u10_col)
  } else {
    if (!"u2" %in% names(df))
      stop("wind = 'present' but there is no `u2` column.")
    out$u2  <- df$u2
    applied <- c(applied, "u2")
  }
  
  # --- unmapped source columns ---------------------------------
  leftover <- setdiff(names(df), c(applied, "date"))
  if (length(leftover)) {
    if (keep_extra) {
      out <- cbind(out, df[, leftover, drop = FALSE])
      message("standardise_forcing(): carried through unmapped: ",
              paste(leftover, collapse = ", "))
    } else {
      message("standardise_forcing(): dropped unmapped: ",
              paste(leftover, collapse = ", "))
    }
  }
  
  validate_forcing(out, strict = strict)
  arrange(out, date)
}


# ---------------------------------------------------------------
# 4. Validation
# ---------------------------------------------------------------

#' Check the canonical schema is satisfied and values are physically plausible
#'
#' The range check exists to catch unit errors, not to judge the climate.
#' A Tmax_C of 305 means the Kelvin conversion did not happen; an
#' ssrd_MJ_m2_day of 300 means running accumulations were summed. Both
#' are silent otherwise and both survive to the results.
validate_forcing <- function(df, strict = FALSE) {
  
  say <- if (strict) stop else warning
  
  req  <- CANONICAL_SCHEMA$name[CANONICAL_SCHEMA$required]
  miss <- setdiff(req, names(df))
  if (length(miss))
    stop("Canonical schema incomplete. Missing: ",
         paste(miss, collapse = ", "))
  
  for (i in seq_len(nrow(CANONICAL_SCHEMA))) {
    nm <- CANONICAL_SCHEMA$name[i]
    if (!nm %in% names(df) || is.na(CANONICAL_SCHEMA$plausible_lo[i])) next
    v   <- df[[nm]]
    lo  <- CANONICAL_SCHEMA$plausible_lo[i]
    hi  <- CANONICAL_SCHEMA$plausible_hi[i]
    bad <- sum(v < lo | v > hi, na.rm = TRUE)
    if (bad > 0)
      say(bad, " value(s) of `", nm, "` outside [", lo, ", ", hi, "] ",
          CANONICAL_SCHEMA$unit[i], ". Observed range: ",
          round(min(v, na.rm = TRUE), 2), " to ",
          round(max(v, na.rm = TRUE), 2),
          ". This usually means a unit conversion was missed.")
  }
  
  # --- internal consistency ------------------------------------
  n <- sum(df$Tmin_C > df$Tmax_C, na.rm = TRUE)
  if (n > 0) say(n, " day(s) with Tmin > Tmax — max/min columns are swapped.")
  
  n <- sum(df$Tdmean_C > (df$Tmax_C + df$Tmin_C) / 2, na.rm = TRUE)
  if (n > 0) say(n, " day(s) with Tdew > Tmean — humidity input is wrong.")
  
  if (any(duplicated(df$date)))
    say(sum(duplicated(df$date)), " duplicated date(s).")
  
  gaps <- as.integer(diff(sort(df$date)))
  if (length(gaps) && any(gaps > 1))
    say(sum(gaps > 1), " gap(s) in the daily series, largest ",
        max(gaps), " days.")
  
  invisible(df)
}