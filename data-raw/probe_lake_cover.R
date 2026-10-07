# ===============================================================
# data-raw/probe_lake_cover.R
#
# Stage C.1, first step: static `lake_cover` probe for the ERA5-Land
# FLake decision (see Ataturk_Master_Checklist_v2.md, Stage C.1 and
# PROJECT_HANDOVER.md section 7).
#
# lake_cover is time-invariant, so this is a small, fast request -
# confirmed against the live CDS "Show API request code" output for
# this exact variable (2 Sep 2026):
#
#   import cdsapi
#   dataset = "reanalysis-era5-land"
#   request = {
#       "variable": ["lake_cover"],
#       "data_format": "netcdf",
#       "download_format": "unarchived"
#   }
#   client = cdsapi.Client()
#   client.retrieve(dataset, request).download()
#
# IMPORTANT — confirmed directly in the CDS request builder: for
# invariant variables the Year/Month/Day/Time AND Geographical area
# controls are all disabled/inert ("Extraction does not apply to time
# invariant variables"). There is no way to request a cropped area for
# lake_cover - the file that comes back is the GLOBAL grid at ERA5-Land's
# native 0.1 degree resolution. It's still small (one 2D field, ~a few
# tens of MB), but crop it locally after download (done below) rather
# than expecting the API to do it.
#
# PREREQUISITE — new CDS platform, personal access token (not the old
# UID:key pair):
#   1. Log in at https://cds.climate.copernicus.eu
#   2. Accept the ERA5-Land licence on the dataset page if you haven't
#      already (Download tab -> Terms of use) - retrieve() fails
#      without this, separately from the API key itself.
#   3. Copy your personal access token from your profile page.
#   4. wf_set_key(key = "<token>") once, interactively - do not commit
#      the token to the repo or paste it into a script.
# ===============================================================

library(ecmwfr)
library(terra)

# --- 1. submit and download ----------------------------------------

request <- list(
  dataset_short_name = "reanalysis-era5-land",
  variable            = "lake_cover",
  data_format         = "netcdf",
  download_format     = "unarchived",
  target              = "lake_cover_global.nc"
)

out_dir <- here::here("data-raw")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

file <- wf_request(
  request = request,
  transfer = TRUE,
  path     = out_dir
)

# --- 2. crop to the Ataturk bounding box and inspect ----------------
# Box from Ataturk_Master_Checklist_v2.md / PROJECT_HANDOVER.md:
# [N, W, S, E] = [38.1, 38.05, 37.35, 39.05]

r   <- rast(file.path(out_dir, "lake_cover_global.nc"))
box <- ext(38.05, 39.05, 37.35, 38.1)   # terra wants (xmin, xmax, ymin, ymax)
r_crop <- crop(r, box)

cat("lake_cover over the Ataturk box:\n")
print(r_crop)
cat("\nValue range:", paste(round(range(values(r_crop), na.rm = TRUE), 4), collapse = " to "), "\n")
cat("Any non-zero cells (FLake active)?", any(values(r_crop) > 0, na.rm = TRUE), "\n")

writeRaster(r_crop, file.path(out_dir, "lake_cover_ataturk_box.tif"), overwrite = TRUE)

# --- Decision point -------------------------------------------------
# cl > 0 somewhere in the box: FLake is running there. lmlt/lmld/ltlt/
#   lblt are live and Stage D can use the FLake-derived route.
# cl == 0 everywhere: FLake is not running on this reservoir (common -
#   FLake's built-in lake mask often misses smaller/narrower reservoirs).
#   Stage D falls back to de Bruin/McJannet equilibrium temperature,
#   with MODIS LST as the primary constraint rather than a validation
#   target. Worth knowing before submitting the multi-day bulk request,
#   since it changes what's worth pulling in that request (still pull
#   the lake_* variables regardless - "cl == 0" doesn't guarantee every
#   cell is inert, and it costs nothing to have them).
