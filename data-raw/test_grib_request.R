# ===============================================================
# data-raw/test_grib_request.R
#
# ONE-OFF DIAGNOSTIC - run this before re-running download_era5land_bulk.R.
#
# Two NetCDF submissions of the bulk request have now failed with the
# same "cost limits exceeded" 403, at very different request sizes
# (105,120 items, then 8,928 items) - so this isn't the item-count
# limit from before. Research (2 Sep 2026, sourced, not recalled)
# turned up a second, separate restriction ECMWF added 3 Apr 2025:
# NetCDF output on ERA5/ERA5-Land now carries its own, much stricter
# and reportedly unquantifiable cost check, because NetCDF requires
# converting the data from its native GRIB format server-side. GRIB
# requests are not subject to this check.
#   https://forum.ecmwf.int/t/limitation-change-on-netcdf-era5-requests/12477
#
# Rather than rewrite all 792 requests in the bulk script a third time
# on an unverified assumption, this script tests the fix on the
# smallest possible case first:
#   - one variable (2m_temperature), one month (Jan 1993), the same box
#   - requests GRIB instead of NetCDF
#   - tries to read the result with terra, since terra/GDAL usually
#     (not confirmed for your GDAL build) reads GRIB directly without
#     a separate conversion step
#
# If this succeeds end-to-end, the bulk script gets the same
# data_format change and is safe to resubmit. If terra can't read the
# GRIB output, we'll need a conversion step (grib_to_netcdf via
# ecCodes, or cdo) before the R pipeline can use these files - better
# to find that out on one small file than after a multi-day 792-request
# queue.
# ===============================================================

library(ecmwfr)
library(terra)

box <- c(38.1, 38.05, 37.35, 39.05)   # N, W, S, E

request <- list(
  dataset_short_name = "reanalysis-era5-land",
  variable            = "2m_temperature",
  year                = "1993",
  month               = "01",
  day                 = sprintf("%02d", 1:31),
  time                = sprintf("%02d:00", 0:23),
  area                = box,
  data_format         = "grib",
  download_format     = "unarchived",
  target              = "test_t2m_199301.grib"
)

out_dir <- here::here("data-raw")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

cat("Submitting GRIB test request (1 variable, 1 month, same box)...\n")

file <- wf_request(
  request  = request,
  transfer = TRUE,
  path     = out_dir
)

cat("\nDownload succeeded:", file, "\n")
cat("Attempting to read it with terra::rast()...\n\n")

r <- rast(file.path(out_dir, "test_t2m_199301.grib"))
print(r)
cat("\nNumber of layers (should be 744, one per hour in a 31-day month):", nlyr(r), "\n")
cat("First-layer value range over the box:\n")
print(range(values(r[[1]]), na.rm = TRUE))

cat("\n=== If the above printed cleanly with 744 layers and sane\n",
    "temperature values (Kelvin, roughly 250-320), GRIB + terra works\n",
    "and the bulk script just needs data_format switched to grib.\n",
    "If rast() errored or the layer count/values look wrong, report\n",
    "the exact output back before we touch the bulk script.\n", sep = "")
