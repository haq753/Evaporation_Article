# ===============================================================
# Daily aggregation and day-of-year climatology
#
# These are two different operations and the old mean_daily() conflated
# them. Keep them apart:
#
#   to_daily()           hourly -> daily time series. One row per date.
#                        This is what the Penman chain consumes.
#
#   daily_climatology()  daily series -> 365-row seasonal cycle.
#                        Display only. Never feed this to a trend test,
#                        a correlation, or a sensitivity analysis.
# ===============================================================

library(dplyr)
library(lubridate)


# ---------------------------------------------------------------
# Hourly -> daily
# ---------------------------------------------------------------

#' Aggregate hourly forcing to a daily time series
#'
#' Variables are aggregated by *type*, not uniformly. Averaging an
#' accumulated flux, or taking the mean of temperature when you need
#' Tmax, are both silent errors that survive all the way to the results.
#'
#' ERA5-Land accumulation convention: accumulated fields reset at 00 UTC
#' and the value stamped 00:00 on day D+1 is the total for day D. If you
#' pass already-deaccumulated hourly rates, set `accum_is_rate = TRUE`
#' and they will be summed instead.
#'
#' @param df    hourly data with a POSIXct `datetime` column
#' @param mean_vars   columns to average (state variables: d2m, sp, u10, v10, lmlt, lmld)
#' @param sum_vars    accumulated fluxes (ssrd, str, strd, tp, e, pev)
#' @param minmax_vars columns needing daily max and min as well as mean (t2m)
#' @param accum_is_rate TRUE if `sum_vars` are already hourly increments
#' @return one row per date
to_daily <- function(df,
                     mean_vars   = character(),
                     sum_vars    = character(),
                     minmax_vars = character(),
                     accum_is_rate = FALSE,
                     tz = "UTC") {
  
  if (!"datetime" %in% names(df))
    stop("`df` needs a POSIXct `datetime` column.")
  
  all_vars <- c(mean_vars, sum_vars, minmax_vars)
  missing  <- setdiff(all_vars, names(df))
  if (length(missing))
    stop("Columns not found: ", paste(missing, collapse = ", "))
  if (length(all_vars) == 0L)
    stop("No variables specified. Nothing to aggregate.")
  
  dup <- all_vars[duplicated(all_vars)]
  if (length(dup))
    stop("Variable(s) listed under more than one aggregation type: ",
         paste(unique(dup), collapse = ", "))
  
  df <- df %>% arrange(datetime)
  
  # --- accumulated fields -------------------------------------------
  # Handled first because they need the D+1 00:00 stamp, not the D dates.
  acc <- NULL
  if (length(sum_vars)) {
    if (accum_is_rate) {
      acc <- df %>%
        mutate(date = as.Date(datetime, tz = tz)) %>%
        group_by(date) %>%
        summarise(across(all_of(sum_vars), ~ sum(.x, na.rm = TRUE)),
                  n_hours = n(), .groups = "drop")
    } else {
      # The running total stamped 00:00 belongs to the PREVIOUS day.
      acc <- df %>%
        filter(hour(datetime) == 0) %>%
        mutate(date = as.Date(datetime, tz = tz) - 1) %>%
        select(date, all_of(sum_vars))
      
      if (nrow(acc) == 0L)
        stop("No 00:00 timestamps found. Either the data is not hourly, ",
             "or it is already deaccumulated - set accum_is_rate = TRUE.")
    }
  }
  
  # --- state variables and min/max -----------------------------------
  base <- df %>%
    mutate(date = as.Date(datetime, tz = tz)) %>%
    group_by(date)
  
  out <- base %>%
    summarise(
      n_hours = n(),
      across(all_of(mean_vars),   ~ mean(.x, na.rm = TRUE)),
      across(all_of(minmax_vars), list(mean = ~ mean(.x, na.rm = TRUE),
                                       max  = ~ max(.x,  na.rm = TRUE),
                                       min  = ~ min(.x,  na.rm = TRUE)),
             .names = "{.col}_{.fn}"),
      .groups = "drop"
    )
  
  if (!is.null(acc)) {
    if (!accum_is_rate) out <- out %>% select(-n_hours)
    out <- left_join(out, acc, by = "date")
  }
  
  # --- completeness flag ---------------------------------------------
  if ("n_hours" %in% names(out)) {
    short <- sum(out$n_hours < 24, na.rm = TRUE)
    if (short > 0)
      warning(short, " day(s) built from fewer than 24 hours. ",
              "Inspect `n_hours` before using these rows.")
  }
  
  arrange(out, date)
}

