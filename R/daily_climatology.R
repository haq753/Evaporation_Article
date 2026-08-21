
# ---------------------------------------------------------------
# Daily -> day-of-year climatology
# ---------------------------------------------------------------

#' Mean seasonal cycle across years
#'
#' FOR PLOTTING ONLY. Collapsing years destroys the interannual
#' variation that every trend and correlation in this study depends on.
#'
#' Leap years are handled by aligning on calendar date rather than raw
#' day-of-year: without this, doy 60 is 29 Feb in a leap year and 1 Mar
#' otherwise, and every date after February is smeared across two
#' neighbouring days of the climatology.
#'
#' @param df   daily data with a `date` column
#' @param vars columns to summarise; default is all numeric columns
#' @param leap "drop" removes 29 Feb and realigns (default);
#'             "raw" uses lubridate::yday() unchanged
#' @return 365 rows (or 366 under "raw"), with n and sd alongside each mean
daily_climatology <- function(df, vars = NULL, leap = c("drop", "raw")) {
  
  leap <- match.arg(leap)
  
  if (!"date" %in% names(df)) stop("`df` needs a `date` column.")
  df$date <- as.Date(df$date)
  
  if (is.null(vars)) {
    vars <- names(df)[vapply(df, is.numeric, logical(1))]
    # Drop bookkeeping columns: a "mean year" is meaningless.
    vars <- setdiff(vars, c("year", "month", "day", "doy", "n_hours"))
    if (length(vars) == 0L)
      stop("No numeric variables left to summarise after dropping ",
           "bookkeeping columns. Pass `vars` explicitly.")
  } else {
    missing <- setdiff(vars, names(df))
    if (length(missing))
      stop("Columns not found: ", paste(missing, collapse = ", "))
  }
  
  if (leap == "drop") {
    n_feb29 <- sum(month(df$date) == 2 & day(df$date) == 29, na.rm = TRUE)
    df <- df %>%
      filter(!(month(date) == 2 & day(date) == 29)) %>%
      mutate(doy = yday(date) - as.integer(leap_year(date) & yday(date) > 59))
    if (n_feb29 > 0)
      message("daily_climatology(): dropped ", n_feb29,
              " leap days and realigned the calendar to 365 days.")
  } else {
    df <- df %>% mutate(doy = yday(date))
  }
  
  out <- df %>%
    group_by(doy) %>%
    summarise(
      n_years = sum(!is.na(.data[[vars[1]]])),
      across(all_of(vars),
             list(mean = ~ mean(.x, na.rm = TRUE),
                  sd   = ~ sd(.x,   na.rm = TRUE)),
             .names = "{.col}_{.fn}"),
      .groups = "drop"
    ) %>%
    arrange(doy)
  
  thin <- sum(out$n_years < 0.8 * max(out$n_years, na.rm = TRUE), na.rm = TRUE)
  if (thin > 0)
    warning(thin, " day(s) of year built from noticeably fewer years ",
            "than the rest. Check `n_years`.")
  
  out
}