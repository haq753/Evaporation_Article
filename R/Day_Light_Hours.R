# ===============================================================
# Daylight hours (FAO-56 Eq. 34)
#
# NOTE ON PLACEMENT
# The two helpers below (solar_declination, sunset_hour_angle) belong
# in R/utils_physics.R, not here. They are reproduced in this file only
# so it runs standalone during the transition. Once utils_physics.R
# exists, delete them from here and from Ra_calc inside ETo_FAO56 and
# Eopen_Penman, so all three functions share one declination.
# ===============================================================


#' Daylight hours, FAO-56 Eq. 34: N = (24 / pi) * omega_s
#'
#' Returns a plain numeric vector, matching the convention of
#' wind10_to_wind2(): the caller decides what to bind it to. This keeps
#' the function agnostic about column names.
#'
#' @param lat_deg  single latitude, decimal degrees, north positive
#' @param dates    vector of Dates, or anything as.Date() accepts
#' @return numeric vector of daylight hours, same length as `dates`
#'
#' @examples
#' # Ataturk, 37.5 N
#' daylight_hours(37.5, as.Date(c("2020-06-21", "2020-12-21")))
#' # equator on the equinox -> ~12.0
#' daylight_hours(0, as.Date("2020-03-21"))
daylight_hours <- function(lat_deg, dates) {
  
  if (length(lat_deg) != 1L || !is.numeric(lat_deg) || is.na(lat_deg))
    stop("lat_deg must be a single non-missing numeric value.")
  if (abs(lat_deg) > 90)
    stop("lat_deg must lie within [-90, 90]. Got ", lat_deg, ".")
  
  dates <- as.Date(dates)
  if (length(dates) == 0L)
    stop("`dates` is empty.")
  if (all(is.na(dates)))
    stop("`dates` contains no parseable dates.")
  
  if (abs(lat_deg) > 66.5)
    warning("lat_deg = ", lat_deg, " is inside the polar circles. ",
            "FAO-56 Eq. 25 is clamped there and returns 0 or 24 h ",
            "rather than a true polar day/night length.")
  
  # Vectorised: no sapply, so no question of the Date class surviving
  # the loop, and no per-element function call overhead.
  J  <- as.integer(strftime(dates, "%j"))
  ws <- sunset_hour_angle(lat_deg, J)
  
  24 / pi * ws
}