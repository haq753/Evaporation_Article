# ===============================================================
# data-raw/wind_station_check.R   (8 Oct 2026)
#
# Question: is ERA5-Land (reservoir-mean 10 m wind, 1.99 m/s) or NASA POWER /
# MERRA-2 (WS10M, 3.53 m/s) closer to OBSERVED wind in this region?
#
# Method: NOAA ISD-Lite hourly station wind (m/s x 10, missing -9999;
# format checked against the NCEI ISD-Lite format document, 8 Oct 2026) for
# stations within RADIUS_KM of the reservoir, reduced to daily means the
# same way as the ERA5-Land series (mean of hourly speed), then compared to
#   (a) POWER WS10M fetched AT THE STATION (like-for-like, same location)
#   (b) the reservoir-mean ERA5-Land series and the reservoir POWER series
#       (regional comparison only: the stations are not on the lake)
#
# Reads: data/daily_forcing_reservoir_1993-2025.csv, data/nasa_power_daily_1993-2025.csv
# Needs network: NCEI (isd-history + isd-lite) and NASA POWER.
# Writes: data/wind_station_summary.csv, data/wind_station_daily.csv
# Cache:  data-raw/isd/   (git-ignored via data-raw/*)
#
# LIMITS, stated up front: ISD does not record anemometer height (WMO
# standard is 10 m, not guaranteed); stations are on land, often airports;
# ERA5-Land at the station is NOT included here (the local ERA5-Land
# download covers only the reservoir box).
# ===============================================================
suppressPackageStartupMessages({ library(dplyr); library(httr); library(nasapower) })

RES_LON   <- 38.5878
RES_LAT   <- 37.6032
RADIUS_KM <- 150
Y0 <- 1993; Y1 <- 2025
MIN_OBS_PER_DAY   <- 8     # obs in the UTC day ...
MIN_BLOCKS        <- 4     # ... spread over all four 6-h blocks
MIN_VALID_DAYS    <- 3000  # drop stations with less than this many valid days
WS_MAX_PLAUSIBLE  <- 40    # m/s; an hourly mean above this invalidates the UTC day
UA <- "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/124.0 Safari/537.36"
dir_isd <- here::here("data-raw", "isd")
dir.create(dir_isd, recursive = TRUE, showWarnings = FALSE)

# ---- helpers ------------------------------------------------------------
get_file <- function(url, dest) {            # returns HTTP status (0 = failed)
  if (file.exists(dest) && file.size(dest) > 0) return(200L)
  r <- tryCatch(GET(url, user_agent(UA), timeout(120), write_disk(dest, overwrite = TRUE)),
                error = function(e) NULL)
  if (is.null(r)) { if (file.exists(dest)) file.remove(dest); return(0L) }
  sc <- status_code(r)
  if (sc != 200L && file.exists(dest)) file.remove(dest)
  sc
}
haversine_km <- function(lon1, lat1, lon2, lat2) {
  r <- pi / 180; a <- sin((lat2 - lat1) * r / 2)^2 +
    cos(lat1 * r) * cos(lat2 * r) * sin((lon2 - lon1) * r / 2)^2
  2 * 6371 * asin(sqrt(a))
}

# ---- 1. candidate stations ---------------------------------------------
hist_f <- file.path(dir_isd, "isd-history.csv")
sc <- get_file("https://www.ncei.noaa.gov/pub/data/noaa/isd-history.csv", hist_f)
if (sc != 200L) stop("Could not download isd-history.csv (HTTP ", sc, "). ",
                     "Download it in a browser to data-raw/isd/ and re-run.")
ish <- read.csv(hist_f, colClasses = "character")
if (ncol(ish) != 11) stop("isd-history.csv has ", ncol(ish), " columns, expected 11.")
names(ish) <- c("USAF", "WBAN", "NAME", "CTRY", "STATE", "ICAO", "LAT", "LON", "ELEV_M", "BEGIN", "END")
ish <- ish |>
  mutate(LAT = suppressWarnings(as.numeric(LAT)), LON = suppressWarnings(as.numeric(LON)),
         ELEV_M = suppressWarnings(as.numeric(ELEV_M)),
         BEGIN = as.integer(BEGIN), END = as.integer(END)) |>
  filter(!is.na(LAT), !is.na(LON), !(LAT == 0 & LON == 0)) |>
  mutate(dist_km = haversine_km(RES_LON, RES_LAT, LON, LAT)) |>
  filter(dist_km <= RADIUS_KM, BEGIN <= 20050101L, END >= 20200101L) |>
  arrange(dist_km)
cat(sprintf("\n%d ISD stations within %d km with records spanning ~2005..2020:\n", nrow(ish), RADIUS_KM))
print(ish |> select(USAF, WBAN, NAME, CTRY, LAT, LON, ELEV_M, BEGIN, END, dist_km) |>
        mutate(dist_km = round(dist_km)), row.names = FALSE)
if (nrow(ish) == 0) stop("No candidate stations; widen RADIUS_KM.")

# ---- 2. download + reduce to daily means -------------------------------
read_isd_lite <- function(f) {
  d <- tryCatch(read.table(gzfile(f), header = FALSE, colClasses = "integer",
                           col.names = c("year", "month", "day", "hour", "temp", "dewp",
                                         "slp", "wdir", "wspd", "sky", "p1", "p6")),
                error = function(e) NULL)
  d
}
daily_station <- function(usaf, wban) {
  out <- list()
  for (y in Y0:Y1) {
    f  <- file.path(dir_isd, sprintf("%s-%s-%d.gz", usaf, wban, y))
    sc <- get_file(sprintf("https://www.ncei.noaa.gov/pub/data/noaa/isd-lite/%d/%s-%s-%d.gz", y, usaf, wban, y), f)
    if (sc == 0L) { message("  network failure for ", usaf, " ", y, " - skipped"); next }
    if (sc != 200L) next                       # 404: station has no file that year
    d <- read_isd_lite(f)
    if (is.null(d) || !nrow(d)) { message("  unreadable: ", basename(f)); next }
    d$ws <- ifelse(d$wspd == -9999L, NA_real_, d$wspd / 10)
    out[[as.character(y)]] <- d
  }
  if (!length(out)) return(NULL)
  d <- bind_rows(out) |>
    mutate(date = as.Date(sprintf("%04d-%02d-%02d", year, month, day)), block = hour %/% 6L)
  # QC: ISD-Lite carries no quality flags. Inspection of Sanliurfa (172700) showed
  # 13 hourly values of 41-87 m/s, all on one day (2018-02-25), against a station
  # mean of 1.5 m/s: instrument/encoding fault, not weather. Rule: any hourly value
  # above WS_MAX_PLAUSIBLE invalidates that whole UTC day. Counts are printed.
  # The scaling itself is fine (typical values are 1-5 m/s, as documented).
  bad_days <- unique(d$date[!is.na(d$ws) & d$ws > WS_MAX_PLAUSIBLE])
  if (length(bad_days))
    cat(sprintf("  QC: dropped %d day(s) with an hourly wind > %d m/s: %s\n", length(bad_days),
                WS_MAX_PLAUSIBLE, paste(format(head(bad_days, 6)), collapse = ", ")))
  d <- d[!d$date %in% bad_days, ]
  cat(sprintf("  median hourly wind %.1f m/s, 99.9th pct %.1f m/s (scaling sanity)\n",
              median(d$ws, na.rm = TRUE), quantile(d$ws, 0.999, na.rm = TRUE)))
  d |>
    group_by(date) |>
    summarise(n_obs = sum(!is.na(ws)), n_blk = n_distinct(block[!is.na(ws)]),
              obs_ws = mean(ws, na.rm = TRUE), .groups = "drop") |>
    filter(n_obs >= MIN_OBS_PER_DAY, n_blk >= MIN_BLOCKS)
}

get_power_ws <- function(usaf, wban, lon, lat) {
  f <- file.path(dir_isd, sprintf("power_ws10m_%s-%s.csv", usaf, wban))
  if (!file.exists(f)) {
    p <- tryCatch(get_power(community = "AG", lonlat = c(lon, lat), pars = "WS10M",
                            dates = c(sprintf("%d-01-01", Y0), sprintf("%d-12-31", Y1)),
                            temporal_api = "daily"),
                  error = function(e) { message("  POWER failed: ", conditionMessage(e)); NULL })
    if (is.null(p)) return(NULL)
    write.csv(data.frame(date = as.Date(p$YYYYMMDD), pow_stn = p$WS10M), f, row.names = FALSE)
    Sys.sleep(1)
  }
  x <- read.csv(f); x$date <- as.Date(x$date); x
}

# ---- 3. regional series already on disk ---------------------------------
era <- read.csv(here::here("data", "daily_forcing_reservoir_1993-2025.csv"))
pwr <- read.csv(here::here("data", "nasa_power_daily_1993-2025.csv"))
reg <- data.frame(date = as.Date(era$date), era_res = era$wspd10_mean_val,
                  pow_res = pwr$WS10M[match(era$date, pwr$date)])

# ---- 4. loop over stations ----------------------------------------------
all_daily <- list(); summ <- list()
for (i in seq_len(nrow(ish))) {
  s <- ish[i, ]
  cat(sprintf("\n[%d/%d] %s-%s %s (%s, %.0f km, %.0f m)\n", i, nrow(ish), s$USAF, s$WBAN, s$NAME, s$CTRY, s$dist_km, s$ELEV_M))
  dd <- daily_station(s$USAF, s$WBAN)
  if (is.null(dd) || nrow(dd) < MIN_VALID_DAYS) {
    cat("  dropped: ", if (is.null(dd)) 0 else nrow(dd), " valid days (< ", MIN_VALID_DAYS, ")\n", sep = ""); next
  }
  pw <- get_power_ws(s$USAF, s$WBAN, s$LON, s$LAT)
  if (is.null(pw)) next
  m <- dd |> inner_join(pw, by = "date") |> inner_join(reg, by = "date") |>
    filter(!is.na(obs_ws), !is.na(pow_stn), !is.na(era_res), !is.na(pow_res))
  m$station <- paste0(s$USAF, "-", s$WBAN); m$name <- s$NAME
  m$month <- as.integer(format(m$date, "%m"))
  m$season <- c("DJF","DJF","MAM","MAM","MAM","JJA","JJA","JJA","SON","SON","SON","DJF")[m$month]
  all_daily[[m$station[1]]] <- m
  summ[[m$station[1]]] <- data.frame(
    station = m$station[1], name = s$NAME, dist_km = round(s$dist_km), elev_m = s$ELEV_M,
    n_days = nrow(m), first = min(m$date), last = max(m$date),
    obs = round(mean(m$obs_ws), 2), power_at_stn = round(mean(m$pow_stn), 2),
    power_res = round(mean(m$pow_res), 2), era5l_res = round(mean(m$era_res), 2),
    ratio_powerStn_obs = round(mean(m$pow_stn) / mean(m$obs_ws), 2),
    ratio_era5lRes_obs = round(mean(m$era_res) / mean(m$obs_ws), 2),
    ratio_powerRes_obs = round(mean(m$pow_res) / mean(m$obs_ws), 2),
    r_powerStn = round(cor(m$obs_ws, m$pow_stn), 2),
    r_era5lRes = round(cor(m$obs_ws, m$era_res), 2),
    r_powerRes = round(cor(m$obs_ws, m$pow_res), 2))
}
if (!length(summ)) stop("No station survived the filters; see messages above.")
S <- bind_rows(summ); D <- bind_rows(all_daily)

cat("\n=== PER-STATION SUMMARY (daily means over days where obs, POWER and ERA5-Land all exist) ===\n")
print(S, row.names = FALSE, width = 200)

cat("\n=== SEASONAL RATIOS vs OBSERVED (mean of series / mean of obs) ===\n")
print(D |> group_by(station, name, season) |>
        summarise(n = n(), obs = round(mean(obs_ws), 2),
                  powerStn_over_obs = round(mean(pow_stn) / mean(obs_ws), 2),
                  era5lRes_over_obs = round(mean(era_res) / mean(obs_ws), 2),
                  powerRes_over_obs = round(mean(pow_res) / mean(obs_ws), 2), .groups = "drop"),
      n = 100, width = 200)

cat("\n=== 5-YEAR MEAN OBSERVED WIND (homogeneity check: instrument/site changes show as steps) ===\n")
print(D |> mutate(bin = 5 * ((as.integer(format(date, "%Y")) - Y0) %/% 5) + Y0) |>
        group_by(station, bin) |> summarise(obs = round(mean(obs_ws), 2), n = n(), .groups = "drop") |>
        tidyr::pivot_wider(names_from = bin, values_from = c(obs), id_cols = station), width = 200)

write.csv(S, here::here("data", "wind_station_summary.csv"), row.names = FALSE)
write.csv(D, here::here("data", "wind_station_daily.csv"), row.names = FALSE)
cat("\nWrote data/wind_station_summary.csv and data/wind_station_daily.csv\n")
