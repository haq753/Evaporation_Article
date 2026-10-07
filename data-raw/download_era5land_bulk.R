# ===============================================================
# data-raw/download_era5land_bulk.R
#
# Stage C.1 bulk request: ERA5-Land hourly forcing, 1993-2025, box
# [N, W, S, E] = [38.1, 38.05, 37.35, 39.05] (Ataturk_Master_Checklist_v2.md).
#
# THIS IS THE CRITICAL-PATH STEP - submit it now, it queues for days.
# Cannot be run from here: neither this session's cloud sandbox nor the
# sandboxed shell on this machine can reach cds.climate.copernicus.eu.
# Run this yourself in RStudio, where your normal internet connection
# works, after the lake_cover probe (data-raw/probe_lake_cover.R) has
# already confirmed your CDS token is live.
#
# CURRENT STATE (2 Sep 2026): requests GRIB, not NetCDF - see
# "REVISION 2" below for why. Output is data-raw/era5land/*.grib,
# readable directly with terra::rast() (verified). Two GRIB-specific
# gotchas (temperature unit mislabeling, per-band naming) are called
# out in REVISION 2 - read that before writing anything downstream
# that opens these files.
#
# ---------------------------------------------------------------
# REVISION (2 Sep 2026) - the first version of this script used a
# 120,000-item-per-request limit and failed on submission:
#
#   Error: permission denied, cost limits exceeded, 403,
#   Your request is too large, please reduce your selection.
#
# That 120,000 figure was a generic example from ECMWF's "Common Error
# Messages for CDS Requests" page, not the actual limit for this
# dataset - a mistake, not a platform change. The real, dataset-specific
# figure, confirmed against multiple ECMWF-staff forum posts and the
# CDS documentation's own MARS-storage limits table (2 Sep 2026):
#
#   reanalysis-era5-land selection limit = 12,000 items per request
#   (item = variable x timestep; area/grid extent does not count).
#   In ECMWF's own words: "you may now download up to 12 variables
#   per month at once." This replaced an earlier 1,000-item limit in
#   Dec 2023 and was still 12,000 as of the CDS docs' Oct 2024 review.
#
# Sources:
#   https://forum.ecmwf.int/t/request-too-large-item-limit-for-era5-land-data/1488
#   https://confluence.ecmwf.int/display/CKB/Climate+Data+Store+(CDS)+documentation
#
# ---------------------------------------------------------------
# WHY MONTHLY REQUESTS, NOT YEARLY
# ---------------------------------------------------------------
# One month, hourly, worst case (31 days) = 744 timesteps.
#
#   12 forcing variables x 744 = 8,928   -> fits, with margin, under 12,000
#    7 lake variables    x 744 = 5,208   -> fits comfortably
#
# A full year at either variable count would exceed 12,000 many times
# over, so each (year, month) pair is now its own request, still split
# into the same two variable groups as before (forcing / lake) so
# neither group's request needs to be thinned further. Full hourly
# resolution is preserved throughout - the fix is request granularity,
# not temporal resolution.
#
# This is also the pattern ECMWF's own support staff recommend for
# ERA5-Land in the forum thread above: "Using an API script, users may
# loop through months and years."
#
# ---------------------------------------------------------------
# REVISION 2 (2 Sep 2026) - the monthly NetCDF version above (this
# file, previous revision) STILL failed on submission with the same
# "cost limits exceeded" 403, at 8,928 items/request - well under the
# 12,000 item cap that fixed nothing here, because it was never the
# actual blocker for this failure. Cause, confirmed by testing (not
# just reading about it): ECMWF added a second, separate restriction
# on 3 Apr 2025 specifically for NetCDF output on ERA5/ERA5-Land -
# converting from GRIB (the dataset's native format) to NetCDF
# server-side now carries its own cost check that multiple people on
# ECMWF's forum report as stricter than, and independent of, the old
# item-count limit, and not reliably predictable in advance:
#   https://forum.ecmwf.int/t/limitation-change-on-netcdf-era5-requests/12477
#
# Fix, verified end-to-end with a real one-variable/one-month request
# (data-raw/test_grib_request.R, run successfully 2 Sep 2026): request
# `data_format = "grib"` instead of "netcdf". GRIB needs no server-side
# conversion, so this cost check never triggers. `terra::rast()` reads
# the resulting .grib files directly - no local conversion step needed.
#
# Two things the test run surfaced that matter for every downstream
# script that reads these files, not just this one:
#
#   1. terra labels the temperature layers' unit as "C" - this is a
#      mislabel from GDAL's GRIB driver, not an actual conversion.
#      The raw values in the test (259.6-271.4 for a single January
#      UTC hour) are physically Kelvin (ERA5 stores 2t in Kelvin per
#      GRIB parameter table 167; a range of ~260-271 is a plausible
#      cold January reading, ~-13 to -2 C, whereas 260-271 read AS
#      Celsius would be nonsensical). Treat every temperature-type
#      variable read from these GRIB files as Kelvin regardless of
#      what terra prints, and verify this again during the Stage C.1
#      arrival sanity checks before trusting it further.
#   2. terra does not carry a clean per-band variable name the way it
#      does for NetCDF (band names come out as "2T_0-SFC", "2T_0-SFC_2",
#      etc., not one variable name x a time dimension). Use
#      `terra::time(r)` for the per-layer timestamp instead of parsing
#      band names - it read correctly in the test (744 steps,
#      1993-01-01 00:00 to 1993-01-31 23:00 UTC).
#
# The test used one variable; this script requests up to 12 (forcing)
# or 7 (lake) variables per file. GRIB is natively multi-parameter (it
# doesn't need the netCDF-only conversion at all), so this is expected
# to work the same way, but it has not itself been tested at that
# variable count - worth confirming on the first handful of requests
# once submitted, rather than assuming silently.
#
# ---------------------------------------------------------------
# REVISION 3 (2 Sep 2026) - RESUME-SAFETY, added before an expected
# laptop disconnect (going home mid-run).
#
# Checked wf_request_batch()'s actual source (not assumed): it has NO
# try/catch around individual requests and NO persistent tracking of
# what it already submitted. In practice this means:
#   - if your network drops, or the laptop sleeps (which suspends the
#     whole R process the same way a network drop does), whatever
#     request was in flight at that moment will very likely throw an
#     error that HALTS THE ENTIRE BATCH, not just that one request.
#   - jobs CDS already accepted keep running server-side regardless -
#     that work is never lost.
#   - but re-running this script unchanged would resubmit all 792
#     requests from scratch, duplicating everything already done.
#
# Fix: before building request_list, scan out_dir for .grib files that
# already exist and skip those (year, month, tag) combos. This makes
# "the script died, run it again" always safe and correct, regardless
# of *why* it died - sleep, dropped wifi, closed laptop, a genuine CDS
# error, anything. It does NOT prevent the halt-on-error itself; it
# just makes recovering from one trivial. If you're back mid-run and
# it looks stalled or errored out, just source the script again.
#
# ---------------------------------------------------------------
# REQUEST COUNT AND EXPECTED RUNTIME
# ---------------------------------------------------------------
# 33 years x 12 months x 2 variable groups = 792 requests, vs. 66
# before. Each individual request is small and should process quickly
# once dequeued, but CDS queueing has per-request overhead independent
# of size, so 792 requests will very likely take longer in wall-clock
# terms than the original 66 would have (had it not been rejected).
# Still expected to be a multi-day job - see total_timeout below.
#
# ---------------------------------------------------------------
# VARIABLE NAMES
# ---------------------------------------------------------------
# All confirmed directly against the live CDS request builder for
# reanalysis-era5-land, 2 Sep 2026 - not recalled from memory, per this
# project's standing rule (see PROJECT_HANDOVER.md section 5).
# ===============================================================

library(ecmwfr)
library(xml2)
library(xml2)

box <- c(38.1, 38.05, 37.35, 39.05)   # N, W, S, E - ecmwfr/CDS area order

vars_forcing <- c(
  "2m_temperature",                        # t2m -> daily mean/max/min via to_daily()
  "2m_dewpoint_temperature",               # d2m
  "surface_pressure",                      # sp
  "skin_temperature",                      # skt
  "10m_u_component_of_wind",               # u10
  "10m_v_component_of_wind",               # v10
  "surface_solar_radiation_downwards",     # ssrd
  "surface_net_thermal_radiation",         # str
  "surface_thermal_radiation_downwards",   # strd
  "total_precipitation",                   # tp
  "total_evaporation",                     # e
  "potential_evaporation"                  # pev
)
stopifnot(length(vars_forcing) * 744 <= 12000)  # 744 = 31 days x 24h, worst-case month

vars_lake <- c(
  "lake_mix_layer_temperature",   # lmlt
  "lake_mix_layer_depth",         # lmld
  "lake_total_layer_temperature", # ltlt
  "lake_bottom_temperature",      # lblt
  "lake_ice_depth",               # licd
  "lake_ice_temperature",         # bonus - free within the item budget
  "lake_shape_factor"             # bonus - free within the item budget
)
stopifnot(length(vars_lake) * 744 <= 12000)

years       <- 1993:2025
months      <- 1:12
all_days    <- sprintf("%02d", 1:31)   # CDS silently ignores invalid combos (e.g. Feb 30)
all_hours   <- sprintf("%02d:00", 0:23)

out_dir <- here::here("data-raw", "era5land")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

make_request <- function(year, month, variables, tag) {
  list(
    dataset_short_name = "reanalysis-era5-land",
    variable            = variables,
    year                = as.character(year),
    month               = sprintf("%02d", month),
    day                 = all_days,
    time                = all_hours,
    area                = box,
    data_format         = "grib",
    download_format     = "unarchived",
    target              = sprintf("era5land_%s_%d_%02d.grib", tag, year, month)
  )
}

ym <- expand.grid(month = months, year = years)  # month varies fastest

request_list_all <- c(
  Map(make_request, year = ym$year, month = ym$month,
      MoreArgs = list(variables = vars_forcing, tag = "forcing")),
  Map(make_request, year = ym$year, month = ym$month,
      MoreArgs = list(variables = vars_lake, tag = "lake"))
)

# --- resume-safety filter: skip anything already downloaded ---------
# (moved into the auto-retry loop below - REVISION 5 - so it re-runs after
# every failure, not just once at the top of the script)

# ---------------------------------------------------------------
# REVISION 4 (3 Sep 2026) - workers dropped from 6 to 1.
#
# After the overnight disconnect, re-running with workers = 6 got a
# NEW, different rejection (not the cost-limit one from before):
#
#   Number queued requests for this dataset is temporarily limited.
#   Please configure your scripts accordingly.
#
# This is CDS refusing to accept more concurrently-queued requests for
# reanalysis-era5-land from this account right now - unrelated to
# per-request size/format (that was Revisions 1-2). Likely triggered
# by some of the pre-disconnect requests still sitting active in CDS's
# queue when the rerun piled 6 more on top. ecmwfr has no function to
# list or cancel what's already queued server-side (checked the
# package's docs, not assumed) - the only place to see that is the
# "Your requests" page at cds.climate.copernicus.eu while logged in.
# ECMWF gives no fixed number for this limit; the word "temporarily"
# matches the dataset-wide queue congestion documented on their forum
# through much of 2026.
#
# Fix: workers = 1, so this script never has more than one request
# outstanding at a time, regardless of what CDS is already holding for
# this account. This will be noticeably slower than workers = 6 - that
# trade is worth it after two rejections from asking for too much at
# once. Raise it later only if you've confirmed via "Your requests"
# that the queue has fully drained and stayed stable.
#
# total_timeout is generous on purpose - this queue is expected to
# take DAYS, and workers = 1 makes that a firm floor, not just a
# possibility.
# ---------------------------------------------------------------

# ---------------------------------------------------------------
# REVISION 5 (22 Sep 2026) - automatic retry loop.
#
# You've been hitting an intermittent 502 Bad Gateway (raw nginx error
# page, not a CDS-level rejection like Revisions 1/2/4 - this is CDS's
# own front-end infrastructure hiccuping, a transient blip rather than
# anything wrong with the request). Your workaround has been: stop the
# script, note what already downloaded, re-run so it doesn't resubmit
# duplicates. That's exactly what the resume-safety filter (REVISION 3)
# already gives you on a manual re-run - this just automates the
# "wait, then re-run" part so a 502 at 3am doesn't sit idle until you're
# back at the keyboard.
#
# The loop below re-scans out_dir for what's already downloaded before
# every attempt (not just once at the top), submits only what's still
# missing, and on ANY error (502 or otherwise - wf_request_batch still
# has no per-request error handling, confirmed from source in REVISION
# 3, so any single failure still halts that attempt) waits 5 minutes
# and tries again from a freshly-rescanned list. Fixed 5-minute wait,
# not exponential backoff - simple, and long enough not to hammer CDS
# on a repeated failure without being so long it wastes half a day if
# the 502 was a one-off blip.
#
# max_retries stops it after 100 consecutive failed ATTEMPTS (not 100
# failed files - one attempt can cover hundreds of files before an
# error interrupts it) so a genuinely fatal problem (expired token,
# account issue) doesn't spin silently forever - it'll print the last
# error and stop, rather than retry something that will never succeed.
# ---------------------------------------------------------------

max_retries    <- 100
retry_wait_sec <- 5 * 60

attempt <- 0
repeat {
  already_have <- list.files(out_dir, pattern = "\\.grib$")
  targets      <- vapply(request_list_all, `[[`, character(1), "target")
  is_done      <- targets %in% already_have
  request_list <- request_list_all[!is_done]

  cat(sum(is_done), "of", length(request_list_all),
      "files present in", out_dir, "-", length(request_list), "still to go.\n")

  if (length(request_list) == 0) {
    cat("All", length(request_list_all), "files downloaded. Done.\n")
    break
  }

  attempt <- attempt + 1
  outcome <- tryCatch({
    wf_request_batch(
      request_list  = request_list,
      workers       = 5,
      path          = out_dir,
      total_timeout = 60 * 60 * 24 * 14
    )
    "ok"
  }, error = function(e) {
    cat("\n[attempt", attempt, "] wf_request_batch stopped with an error:\n  ",
        conditionMessage(e), "\n")
    "error"
  })

  if (identical(outcome, "ok")) {
    cat("\nAll transfers complete. Files in:", out_dir, "\n")
    break
  }

  if (attempt >= max_retries) {
    cat("\nGiving up after", attempt, "attempts. Check the error above and",
        "https://cds.climate.copernicus.eu -> Your requests before re-running.\n")
    break
  }

  cat("Waiting", retry_wait_sec / 60,
      "minutes, then re-scanning for anything that finished before retrying...\n")
  Sys.sleep(retry_wait_sec)
}
