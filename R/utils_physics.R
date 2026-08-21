# ===============================================================
# R/utils_physics.R
#
# Single source of truth for the physics shared by ETo_FAO56() and
# Eopen_Penman(). Nothing in this file knows about column names except
# prepare_forcing(), which is the one place the naming convention is
# enforced.
#
# Source this BEFORE the model functions.
#
# CONVENTIONS PRESERVED FROM FAO-56 — do not "tidy" these:
#   - Eq. 39 uses T + 273.16; the ETo denominator uses T + 273.
#     They differ. Both are reproduced exactly as published.
#   - Tetens coefficients 17.27 / 237.3 are used EVERYWHERE, including
#     for relative humidity. The old Relative_Humidity.R used Magnus
#     (17.625 / 243.04) instead; that file is superseded by rh_from_td()
#     below and should be archived.
# ===============================================================

# --- constants -------------------------------------------------

SIGMA_MJ    <- 4.903e-9   # Stefan-Boltzmann, MJ K^-4 m^-2 d^-1
LAMBDA_MJ   <- 2.45       # latent heat of vaporisation, MJ kg^-1.
# Fixed, not temperature-dependent: the 6.43
# wind constant is expressed in MJ m^-2 d^-1
# kPa^-1 and was derived on this basis.
EPS_WATER   <- 0.97       # water surface emissivity
ALBEDO_WATER <- 0.08
ALBEDO_GRASS <- 0.23


# --- vapour pressure -------------------------------------------

#' Saturation vapour pressure (kPa) at temperature Tc (deg C). FAO-56 Eq. 11.
es_kPa <- function(Tc) 0.6108 * exp(17.27 * Tc / (Tc + 237.3))

#' Slope of the saturation vapour pressure curve (kPa/deg C). FAO-56 Eq. 13.
slope_es <- function(Tc) 4098 * es_kPa(Tc) / (Tc + 237.3)^2

#' Psychrometric constant (kPa/deg C). FAO-56 Eq. 8.
gamma_psy <- function(P_kPa) 0.000665 * P_kPa

#' Actual vapour pressure (kPa) from relative humidity.
ea_from_rh <- function(Tc, RH_pct) es_kPa(Tc) * RH_pct / 100

#' Relative humidity (%) from air and dewpoint temperature.
#' Replaces the Magnus-based compute_RH() in Relative_Humidity.R.
rh_from_td <- function(Tc, Td_C) 100 * es_kPa(Td_C) / es_kPa(Tc)

#' Dewpoint (deg C) from actual vapour pressure. Inverse of es_kPa.
td_from_ea <- function(ea_kPa) {
  if (any(ea_kPa <= 0, na.rm = TRUE))
    stop("ea must be strictly positive.")
  l <- log(ea_kPa / 0.6108)
  237.3 * l / (17.27 - l)
}


# --- atmospheric pressure --------------------------------------

#' Atmospheric pressure (kPa) from elevation. FAO-56 Eq. 7.
#' Prefer measured surface pressure where available.
pressure_from_elev <- function(elev_m) {
  101.3 * ((293 - 0.0065 * elev_m) / 293)^5.26
}


# --- solar geometry --------------------------------------------

#' Solar declination (radians). FAO-56 Eq. 24.
solar_declination <- function(J) 0.409 * sin(2 * pi * J / 365 - 1.39)

#' Inverse relative Earth-Sun distance. FAO-56 Eq. 23.
inverse_earth_sun <- function(J) 1 + 0.033 * cos(2 * pi * J / 365)

#' Sunset hour angle (radians). FAO-56 Eq. 25.
#' Clamped to the acos domain: a guard against misuse poleward of
#' 66.5 deg, not a correction at mid-latitude.
sunset_hour_angle <- function(lat_deg, J) {
  phi   <- lat_deg * pi / 180
  delta <- solar_declination(J)
  acos(pmin(1, pmax(-1, -tan(phi) * tan(delta))))
}

#' Extraterrestrial radiation (MJ m^-2 d^-1). FAO-56 Eq. 21.
Ra_calc <- function(lat_deg, J) {
  lat   <- lat_deg * pi / 180
  dr    <- inverse_earth_sun(J)
  delta <- solar_declination(J)
  ws    <- sunset_hour_angle(lat_deg, J)
  (24 * 60 / pi) * 0.0820 * dr *
    (ws * sin(lat) * sin(delta) + cos(lat) * cos(delta) * sin(ws))
}

#' Clear-sky solar radiation (MJ m^-2 d^-1). FAO-56 Eq. 37.
Rso_calc <- function(Ra, elev_m) (0.75 + 2e-5 * elev_m) * Ra

#' Daylight hours. FAO-56 Eq. 34.
daylight_hours <- function(lat_deg, J) 24 / pi * sunset_hour_angle(lat_deg, J)


# --- net longwave ----------------------------------------------

#' Net longwave using air temperature as a surrogate for surface
#' temperature. FAO-56 Eq. 39.
#'
#' DELIVERABLE NOTE: the clamps below (Rs/Rso floored at 0.25 and capped
#' at 1.0, f_cloud floored at 0.05) are a deliberate deviation from
#' FAO-56 as published. They prevent unphysical Rnl on days with
#' anomalous Rs/Rso. This must be declared in Methods.
Rnl_fao56 <- function(Tmax_C, Tmin_C, ea, Rs, Rso) {
  sigT4   <- SIGMA_MJ * ((Tmax_C + 273.16)^4 + (Tmin_C + 273.16)^4) / 2
  rel_sw  <- pmin(1, pmax(0.25, ifelse(Rso > 0, Rs / Rso, 0.25)))
  f_cloud <- pmax(0.05, 1.35 * rel_sw - 0.35)
  sigT4 * (0.34 - 0.14 * sqrt(pmax(ea, 0))) * f_cloud
}

#' Net longwave for open water from measured/modelled water surface
#' temperature and downwelling longwave.
#'   Rnl = eps_w * (sigma * Tw^4 - Ld)
#' No cloud factor: Ld already carries the cloud signal.
Rnl_water <- function(Tw_C, Ld_MJ, eps_w = EPS_WATER) {
  eps_w * (SIGMA_MJ * (Tw_C + 273.16)^4 - Ld_MJ)
}


# --- shared pre-processing -------------------------------------

REQUIRED_FORCING <- c("date", "Tmax_C", "Tmin_C", "Tdmean_C",
                      "u2", "ssrd_MJ_m2_day")

#' Everything ETo_FAO56() and Eopen_Penman() compute identically.
#'
#' Both models call this, so they cannot drift apart in their treatment
#' of vapour pressure, pressure, or radiation geometry. If a diagnostic
#' ever needs to change, it changes here once.
#'
#' @param df date-sorted forcing data
#' @return list of derived quantities, all same length as nrow(df)
prepare_forcing <- function(df, lat_deg, elev_m, verbose = TRUE) {
  
  miss <- setdiff(REQUIRED_FORCING, names(df))
  if (length(miss))
    stop("Missing required columns: ", paste(miss, collapse = ", "))
  
  J <- as.integer(strftime(df$date, "%j"))
  
  P_kPa <- if ("sp_Pa" %in% names(df)) df$sp_Pa / 1000 else
    pressure_from_elev(elev_m)
  
  Tmean <- (df$Tmax_C + df$Tmin_C) / 2
  
  # FAO-56 Eq. 12: mean of es(Tmax) and es(Tmin), NOT es(Tmean)
  es  <- (es_kPa(df$Tmax_C) + es_kPa(df$Tmin_C)) / 2
  ea  <- es_kPa(df$Tdmean_C)
  vpd <- pmax(0, es - ea)
  
  n_inv <- sum(df$Tdmean_C > Tmean, na.rm = TRUE)
  if (n_inv > 0 && verbose)
    warning(n_inv, " day(s) with Tdew > Tmean - check the humidity input.")
  
  Rs  <- pmax(df$ssrd_MJ_m2_day, 0)
  Ra  <- Ra_calc(lat_deg, J)
  Rso <- Rso_calc(Ra, elev_m)
  
  n_bad <- sum(Rs > Ra, na.rm = TRUE)
  if (n_bad > 0 && verbose)
    warning(n_bad, " day(s) with Rs > Ra - radiation units or ",
            "aggregation are wrong.")
  
  list(J = J, P_kPa = P_kPa, Tmean = Tmean,
       es = es, ea = ea, vpd = vpd,
       Delta = slope_es(Tmean), gamma = gamma_psy(P_kPa),
       Rs = Rs, Ra = Ra, Rso = Rso)
}


