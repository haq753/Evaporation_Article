# ===============================================================
# data-raw/inspect_grib_structure.R
#
# Run this FIRST, before writing/trusting any daily-aggregation code
# against the 792 downloaded GRIB files. It answers a question the
# single-variable test back on 2 Sep 2026 (test_grib_request.R) could
# not: when a GRIB file holds 12 DIFFERENT variables x 744 hours in one
# file, how does terra lay out and name the resulting layers?
#
# Two things are genuinely unverified right now, and both are
# foundational to every number the paper will report:
#
#   1. LAYER ORDER - are the ~8,928 layers in a forcing file grouped
#      by variable (all of t2m's 744 hours, then all of d2m's, ...) or
#      interleaved by time (t2m@h1, d2m@h1, ..., pev@h1, t2m@h2, ...)?
#      Get this wrong and a "t2m" column could silently be full of
#      dewpoint or pressure values instead.
#   2. NAMING - does each of the 12 variables get its own distinct
#      varname (e.g. "2T_0-SFC" for temperature, "2D_0-SFC" for
#      dewpoint, "SP_0-SFC" for pressure...), or does terra's GRIB
#      driver do something less clean with 12 different parameters
#      than it did with the one-variable test?
#
# The single-variable test confirmed terra reads GRIB natively and
# that `terra::time()` gives correct per-layer timestamps - but it
# used ONE variable, so it could not distinguish "grouped by variable"
# from "interleaved by time", and it could not show what happens to
# variable-to-name mapping with 12 variables instead of 1.
#
# Rather than guess and build 200 lines of aggregation logic on top of
# a guess (this project has been burned enough times this month by
# confident-but-wrong assumptions about CDS/ERA5-Land internals to
# know better), this script just prints the real structure so the next
# script can be written against what's actually there.
#
# Run it, then send back the full console output (or at minimum: the
# unique varnames found, how many layers per unique varname, and the
# first 3 and last 3 timestamps for the first two distinct varnames).
# ===============================================================

library(terra)

era_dir <- here::here("data-raw", "era5land")

inspect_one <- function(path) {
  cat("\n=====================================================\n")
  cat("File:", basename(path), "\n")
  r <- rast(path)
  cat("Total layers:", nlyr(r), "\n")
  
  nm <- varnames(r)
  base_nm <- sub("_[0-9]+$", "", nm)          # strip terra's "_2", "_3", ... suffix
  tab <- table(base_nm)
  cat("\nUnique base varnames found (", length(tab), "distinct):\n")
  print(tab)
  
  cat("\n--- Per-varname detail (first 3 / last 3 layer timestamps) ---\n")
  tm <- time(r)
  for (v in names(tab)) {
    idx <- which(base_nm == v)
    cat(sprintf("\n  %-14s n_layers=%-5d first_idx=%-6d last_idx=%d\n",
                v, length(idx), min(idx), max(idx)))
    show_idx <- idx[c(1:min(3, length(idx)))]
    cat("    first timestamps: ", paste(as.character(tm[show_idx]), collapse = " | "), "\n")
    show_idx2 <- idx[max(1, length(idx) - 2):length(idx)]
    cat("    last timestamps:  ", paste(as.character(tm[show_idx2]), collapse = " | "), "\n")
    # is this variable's layer set contiguous (grouped) or spread out (interleaved)?
    contiguous <- all(diff(idx) == 1)
    cat("    layer indices contiguous (grouped by variable)?", contiguous, "\n")
  }
  invisible(r)
}

cat("###################################################\n")
cat("# FORCING FILE (12 variables) - January 1993\n")
cat("###################################################\n")
inspect_one(file.path(era_dir, "era5land_forcing_1993_01.grib"))

cat("\n\n###################################################\n")
cat("# LAKE FILE (7 variables) - January 1993\n")
cat("###################################################\n")
inspect_one(file.path(era_dir, "era5land_lake_1993_01.grib"))

cat("\n\nDone. Paste the full output back - especially the varname table,\n",
    "the 'contiguous' answers, and the first/last timestamps per variable.\n", sep = "")