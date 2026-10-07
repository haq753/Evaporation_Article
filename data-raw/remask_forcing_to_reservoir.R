# ===============================================================
# data-raw/remask_forcing_to_reservoir.R
#
# Stage C.2 (Ataturk_Master_Checklist_v2.md): re-mask the ERA5-Land
# cells to the real reservoir polygon and compute reservoir-mean daily
# forcing + the between-cell spread - the step flagged to "also
# resolve/clarify the licd and low-u2 findings" from Stage C.1.
#
# WHY THIS MATTERS, QUANTIFIED (1 Oct 2026): the request box is 8x10 =
# 80 ERA5-Land cells (0.1 deg each), but the real reservoir (per the
# final JRC polygon) only touches 34 of them, and barely for most of
# those - the sum of per-cell overlap fractions is ~7.67, meaning the
# box is equivalent to fewer than 8 "fully water" cells out of 80. So
# build_daily_forcing.R's box-mean (`global(r, "mean")` over ALL 80
# cells, no masking) is overwhelmingly a LAND-dominated average, not a
# reservoir one - that's the direct, concrete explanation for the u2
# (wind, pulled low by sheltered land cells) and licd (lake ice depth,
# diluted by land/non-FLake cells) findings flagged `[~]` in Stage C.1,
# not just a "working hypothesis" anymore.
#
# TWO DIFFERENT WEIGHTS, DELIBERATELY NOT ONE:
#   - W_JRC: exact fractional overlap of each of the 80 cells with the
#     REAL reservoir extent (the final JRC polygon, original box,
#     confirmed with Aziz 1 Oct 2026). Used for every ordinary
#     atmospheric forcing variable (t2m, d2m, skt, sp, u10/v10/wspd10,
#     ssrd, str, strd, tp, e, pev) - these are physically meaningful in
#     every cell regardless of ECMWF's own lake mask, so weighting by
#     the TRUE water-surface fraction is what "reservoir-mean forcing"
#     should mean.
#   - W_CL: ERA5-Land's OWN internal lake_cover fraction per cell
#     (data-raw/lake_cover_ataturk_box.tif, from probe_lake_cover.R).
#     Used ONLY for the 7 FLake lake-state variables (lmlt, lmld,
#     lblt, ltlt, lshf, lict, licd). These are FLake MODEL outputs,
#     only physically meaningful where FLake itself is actually
#     running (cl > 0) - ECMWF's internal lake mask and the real-world
#     JRC extent don't have to agree (and the Stage C.1 notes already
#     show they don't fully), and weighting FLake's own state variables
#     by the real-world extent instead of by where the model is even
#     active would average in dummy/land-scheme values from cells
#     where FLake isn't running at all.
#
# HOW W_JRC WAS COMPUTED (and why it's a hardcoded literal, not
# computed live by this script): this session hit two real doc-vs-
# reality mismatches already today (DAHITI's documented response
# schema not matching its live output; GRIB variable labels needing a
# verified lookup table rather than assumed names) - rather than trust
# terra::extract(exact=TRUE) or rasterize(cover=TRUE)'s exact output
# format/behavior sight unseen for something this foundational,
# W_JRC was computed independently in Python (shapely), checked, and
# is embedded below as a literal. Method, exactly:
#   1. Read data/ataturk_jrc_max_extent_polygon.geojson (the FINAL,
#      original-box polygon, confirmed correct with Aziz 1 Oct 2026 -
#      NOT the widened-box one superseded earlier that day).
#   2. shapely.validation.make_valid() to fix one tiny ring self-
#      intersection in the exported polygon (confirmed negligible:
#      ~9e-16 deg^2 area difference, i.e. floating-point noise, not a
#      real geometry problem).
#   3. For each of the 80 cells (built directly from the box, 0.1 deg
#      resolution, row 0 = north/top, col 0 = west, matching terra's
#      default row-major cell order): exact polygon-cell intersection
#      area / cell area.
# Regenerate the same way if the polygon is ever rebuilt - don't hand-
# edit these numbers.
#
# GRID CONFIRMED BY DIRECT INSPECTION, NOT ASSUMED FROM THE REQUESTED
# BOX: read era5land_forcing_1993_01.grib, era5land_lake_1993_01.grib,
# and the same pair for 2010-06 and 2025-12 directly with GDAL/rasterio
# (Python) on 1 Oct 2026 - all five files (same grid, four different
# year/months, both forcing and lake) share IDENTICAL geometry: 8 rows
# x 10 cols, 0.1 deg resolution, bounds lon [38.05, 39.05], lat [37.35,
# 38.15]. Note the NORTH edge is 38.15, not the box spec's nominal
# 38.1 - CDS (like terra's own crop()) never splits a cell, so the
# actual returned grid is very slightly larger than the nominal box on
# that side. This matches lake_cover_ataturk_box.tif's own grid
# exactly too. load_grib_reservoir() below checks every one of the 792
# files against this confirmed grid and stops loudly on any mismatch -
# W_JRC's cell order is meaningless against a different grid.
#
# TEST ON 1993 FIRST - RUN_YEARS below is deliberately narrow. Widen to
# 1993:2025 only after checking this year's output (see the comparison
# block at the bottom, which checks against the existing box-mean file
# automatically) - same discipline build_daily_forcing.R itself used.
# ===============================================================

library(terra)
library(dplyr)
library(tidyr)
library(purrr)
library(lubridate)

era_dir <- here::here("data-raw", "era5land")
out_dir <- here::here("data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- verified variable lookup table (base_label/code/kind/units copied
# VERBATIM from build_daily_forcing.R's already-verified VAR_TABLE -
# only the new weight_set column is added) ------------------------
VAR_TABLE <- tribble(
  ~base_label,                                  ~code,   ~kind,           ~unit_in, ~unit_out, ~weight_set,
  "2T_0-SFC",                                    "t2m",   "instantaneous", "K",      "C",       "jrc",
  "2D_0-SFC",                                    "d2m",   "instantaneous", "K",      "C",       "jrc",
  "SKT_0-SFC",                                   "skt",   "instantaneous", "K",      "C",       "jrc",
  "SP_0-SFC",                                    "sp",    "instantaneous", "Pa",     "kPa",     "jrc",
  "10U_0-SFC",                                   "u10",   "instantaneous", "m/s",    "m/s",     "jrc",
  "10V_0-SFC",                                   "v10",   "instantaneous", "m/s",    "m/s",     "jrc",
  "SSRD_0-SFC",                                  "ssrd",  "accumulated",   "J/m2",   "MJ/m2",   "jrc",
  "STR_0-SFC",                                   "str",   "accumulated",   "J/m2",   "MJ/m2",   "jrc",
  "STRD_0-SFC",                                  "strd",  "accumulated",   "J/m2",   "MJ/m2",   "jrc",
  "TP_0-SFC",                                    "tp",    "accumulated",   "m",      "mm",      "jrc",
  "E_0-SFC",                                     "e",     "accumulated",   "m",      "mm",      "jrc",
  "var251 of table 228 of center ECMWF_0-SFC",   "pev",   "accumulated",   "m",      "mm",      "jrc",
  "var8 of table 228 of center ECMWF_0-SFC",     "lmlt",  "instantaneous", "K",      "C",       "cl",
  "var9 of table 228 of center ECMWF_0-SFC",     "lmld",  "instantaneous", "m",      "m",       "cl",
  "var10 of table 228 of center ECMWF_0-SFC",    "lblt",  "instantaneous", "K",      "C",       "cl",
  "var11 of table 228 of center ECMWF_0-SFC",    "ltlt",  "instantaneous", "K",      "C",       "cl",
  "var12 of table 228 of center ECMWF_0-SFC",    "lshf",  "instantaneous", "-",      "-",       "cl",
  "var13 of table 228 of center ECMWF_0-SFC",    "lict",  "instantaneous", "K",      "C",       "cl",
  "var14 of table 228 of center ECMWF_0-SFC",    "licd",  "instantaneous", "m",      "m",       "cl",
)

strip_suffix <- function(x) sub("_[0-9]+$", "", x)

# --- confirmed grid (direct GDAL/rasterio inspection, 1 Oct 2026) ---
GRID_NROW <- 8L; GRID_NCOL <- 10L
GRID_XMIN <- 38.05; GRID_XMAX <- 39.05
GRID_YMIN <- 37.35; GRID_YMAX <- 38.15

check_grid <- function(r, path) {
  e <- ext(r)
  ok <- nrow(r) == GRID_NROW && ncol(r) == GRID_NCOL &&
    abs(e$xmin - GRID_XMIN) < 1e-4 && abs(e$xmax - GRID_XMAX) < 1e-4 &&
    abs(e$ymin - GRID_YMIN) < 1e-4 && abs(e$ymax - GRID_YMAX) < 1e-4
  if (!ok) {
    stop("Grid mismatch in ", basename(path), ": got ", nrow(r), "x", ncol(r),
         " extent (", round(e$xmin,5), ",", round(e$xmax,5), ",",
         round(e$ymin,5), ",", round(e$ymax,5), ") - expected ",
         GRID_NROW, "x", GRID_NCOL, " extent (", GRID_XMIN, ",", GRID_XMAX, ",",
         GRID_YMIN, ",", GRID_YMAX, "). W_JRC/W_CL's cell order would be WRONG ",
         "for this file - STOP, do not let this run silently against a different grid.")
  }
}

# --- W_JRC: exact fractional reservoir overlap per cell, hardcoded ---
# (see header for the full derivation). Row-major from north-west:
# row 0 = northernmost, col 0 = westernmost - matches terra's default
# values() cell order for a raster with ymax at the top.
W_JRC <- c(
  0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000,
  0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.003420,
  0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.013055, 0.000000, 0.000000, 0.000000, 0.104895,
  0.000000, 0.000000, 0.000000, 0.000000, 0.000000, 0.033897, 0.122228, 0.000691, 0.368641, 0.240916,
  0.000000, 0.007240, 0.160798, 0.072521, 0.087679, 0.008234, 0.240624, 0.285725, 0.412432, 0.193312,
  0.000000, 0.087170, 0.494231, 0.625175, 0.201553, 0.605168, 0.297578, 0.232032, 0.029318, 0.006859,
  0.000000, 0.000000, 0.321327, 0.626789, 0.848540, 0.572561, 0.000431, 0.000000, 0.000000, 0.000000,
  0.000000, 0.000000, 0.000000, 0.001928, 0.197013, 0.167927, 0.000000, 0.000000, 0.000000, 0.000000
)
stopifnot(length(W_JRC) == GRID_NROW * GRID_NCOL)

# --- W_CL: ERA5-Land's own lake_cover fraction, read from file (not
# hardcoded) so it can't drift out of sync with probe_lake_cover.R ---
lake_cover_path <- here::here("data-raw", "lake_cover_ataturk_box.tif")
r_cl <- rast(lake_cover_path)
check_grid(r_cl, lake_cover_path)
W_CL <- as.numeric(values(r_cl))

# --- startup diagnostic: print both weight grids reshaped to 8x10 so
# the cell-order assumption above is visible and eyeball-checkable,
# not just trusted from terra's documented (but today, twice-burned-on)
# conventions. Compare this printout's W_JRC grid against the numbers
# in the header comment above - they should match exactly.
cat("=== W_JRC (real reservoir fraction per cell, row 0 = north) ===\n")
print(round(matrix(W_JRC, nrow = GRID_NROW, ncol = GRID_NCOL, byrow = TRUE), 3))
cat("\n=== W_CL (ERA5-Land's own lake_cover fraction per cell, row 0 = north) ===\n")
print(round(matrix(W_CL, nrow = GRID_NROW, ncol = GRID_NCOL, byrow = TRUE), 3))
cat("\nIf these two grids don't show roughly the same general shape/location\n")
cat("(even though the exact values legitimately differ - that's the whole\n")
cat("point of using two different weight sets), STOP and check the cell-order\n")
cat("assumption before trusting anything below.\n\n")

# --- weighted mean + weighted sd across cells, TRUE matrix vectorization
# (2 Oct 2026 rewrite): the original version here was commented
# "vectorized over layers" but wasn't - apply()/vapply() are still
# R-level loops, one iteration per hourly column, and with ~744
# columns x 19 variables x 396 months (~5.6M iterations, each paying
# R's closure-call overhead) this was the dominant reason the full
# 1993-2025 run was taking hours rather than minutes (flagged by Aziz,
# 2 Oct 2026, after the version-mismatch crash was fixed and the run
# had been going ~2h). Confirmed empirically from the 1993 test output
# that NA never occurs in this data (every `*_cellsd` column reported
# "NA count: 0"), but NA-handling is kept anyway via a zeroed-weight
# trick so this isn't silently wrong if a future month ever does have
# a gap - verify, don't assume it stays NA-free forever.
# vmat: ncell(80) x nlayer matrix of raw per-cell values. w: length-80
# weight vector (W_JRC or W_CL). Returns per-layer weighted mean (the
# "reservoir-mean" value) and weighted sd (the between-cell spread).
weighted_cell_stats <- function(vmat, w) {
  sel <- w > 0
  if (!any(sel)) stop("No cells with positive weight - check W_JRC/W_CL.")
  vsel <- vmat[sel, , drop = FALSE]        # ncell_sel x ncol
  wsel <- w[sel]                            # ncell_sel

  na_mask <- is.na(vsel)
  has_na <- any(na_mask)
  vsel0 <- vsel
  if (has_na) vsel0[na_mask] <- 0

  # weight matrix, same shape as vsel0, with NA cells' weight zeroed
  # per-column (so a cell's NA on one hour doesn't poison other hours)
  wmat <- matrix(wsel, nrow = length(wsel), ncol = ncol(vsel))
  if (has_na) wmat[na_mask] <- 0

  wsum <- colSums(wmat)
  wmean <- colSums(vsel0 * wmat) / wsum
  wmean[wsum == 0] <- NA_real_

  meanmat <- matrix(wmean, nrow = nrow(vsel), ncol = ncol(vsel), byrow = TRUE)
  sqdiff <- (vsel0 - meanmat)^2
  if (has_na) sqdiff[na_mask] <- 0
  n_valid <- if (has_na) colSums(!na_mask) else rep(nrow(vsel), ncol(vsel))
  wvar <- colSums(wmat * sqdiff) / wsum
  wsd <- sqrt(wvar)
  wsd[n_valid < 2 | wsum == 0] <- NA_real_

  list(mean = as.numeric(wmean), sd = as.numeric(wsd))
}

# --- load one (year, month, tag) GRIB file, reservoir/cl-weighted ---
# Mirrors build_daily_forcing.R's load_grib_long() for parsing, unit
# conversion and layer-count sanity checks (copied, not shared via
# source(), to avoid touching the already-verified box-mean pipeline).
# The ONLY structural change: spatial reduction is a weighted mean+sd
# over the appropriate cell subset, not global(r, "mean") over all 80.
load_grib_reservoir <- function(year, month, tag) {
  path <- file.path(era_dir, sprintf("era5land_%s_%d_%02d.grib", tag, year, month))
  if (!file.exists(path)) stop("Missing file: ", path)

  r <- rast(path)
  check_grid(r, path)
  base_lab <- strip_suffix(varnames(r))
  tm <- time(r)
  vmat_all <- values(r)   # 80 x nlyr, row-major cell order matching W_JRC/W_CL

  unmatched <- setdiff(unique(base_lab), VAR_TABLE$base_label)
  if (length(unmatched) > 0) {
    stop("unrecognized variable label(s) in ", basename(path), ": ",
         paste(unmatched, collapse = ", "), " - update VAR_TABLE, do not guess.")
  }

  n_hours_expected <- length(unique(tm))
  rows <- list()
  vmat_u10 <- vmat_v10 <- tm_u10 <- tm_v10 <- NULL

  for (i in seq_len(nrow(VAR_TABLE))) {
    lab <- VAR_TABLE$base_label[i]
    idx <- which(base_lab == lab)
    if (length(idx) == 0) next
    code <- VAR_TABLE$code[i]

    vmat <- vmat_all[, idx, drop = FALSE]
    vmat <- switch(paste(VAR_TABLE$unit_in[i], VAR_TABLE$unit_out[i]),
      "K C"        = vmat - 273.15,
      "Pa kPa"     = vmat / 1000,
      "J/m2 MJ/m2" = vmat / 1e6,
      "m mm"       = vmat * 1000,
      vmat)

    if (ncol(vmat) != n_hours_expected) {
      stop("layer-count mismatch in ", basename(path), " for ", code, " - expected ",
           n_hours_expected, " hourly layers, got ", ncol(vmat))
    }

    w <- if (VAR_TABLE$weight_set[i] == "cl") W_CL else W_JRC
    st <- weighted_cell_stats(vmat, w)
    rows[[code]] <- tibble(code = code, kind = VAR_TABLE$kind[i], time = tm[idx],
                            value = st$mean, cellsd = st$sd)

    # Stash the RAW (unconverted, but u10/v10 need no conversion anyway)
    # per-cell matrices for u10/v10 so wspd10 can be derived per-cell
    # below, before any reduction to a scalar - see the header note on
    # why this has to happen at the per-cell level, not from already-
    # reduced u10/v10 means (same reasoning as build_daily_forcing.R's
    # own component-vs-magnitude fix, applied spatially here too).
    if (code == "u10") { vmat_u10 <- vmat; tm_u10 <- tm[idx] }
    if (code == "v10") { vmat_v10 <- vmat; tm_v10 <- tm[idx] }
  }

  if (!is.null(vmat_u10) && !is.null(vmat_v10)) {
    stopifnot(identical(tm_u10, tm_v10))  # same file - verify, don't assume
    vmat_wspd10 <- sqrt(vmat_u10^2 + vmat_v10^2)
    st_w <- weighted_cell_stats(vmat_wspd10, W_JRC)
    rows[["wspd10"]] <- tibble(code = "wspd10", kind = "instantaneous", time = tm_u10,
                                value = st_w$mean, cellsd = st_w$sd)
  }

  bind_rows(rows)
}

# --- accumulated-variable end-of-day extraction (same logic/handling
# of the documented end-of-record gap as build_daily_forcing.R) -----
get_next_day_00h_reservoir <- function(year, month) {
  ny <- year; nm <- month + 1
  if (nm > 12) { ny <- year + 1; nm <- 1 }
  path <- file.path(era_dir, sprintf("era5land_forcing_%d_%02d.grib", ny, nm))
  if (!file.exists(path)) return(NULL)
  df <- load_grib_reservoir(ny, nm, "forcing") |> filter(kind == "accumulated")
  df |> filter(time == min(time)) |> select(code, value, cellsd) |>
    rename(next00 = value, next00_cellsd = cellsd)
}

build_month_reservoir <- function(year, month) {
  # Progress line (2 Oct 2026 addition): with no per-month output the
  # original version gave no way to tell "still working" from "stuck"
  # during a multi-hour run - Aziz hit exactly that after the full
  # 1993-2025 run had been going ~2h with a silent console. cat() (not
  # message()) so it shows up directly in the console/log either way.
  cat(sprintf("[%s] Processing %d-%02d...\n", format(Sys.time(), "%H:%M:%S"), year, month))
  forcing <- load_grib_reservoir(year, month, "forcing")  # wspd10 already included
  lake    <- load_grib_reservoir(year, month, "lake")

  inst <- bind_rows(forcing, lake) |>
    filter(kind == "instantaneous") |>
    mutate(date = as.Date(time)) |>
    group_by(date, code) |>
    summarise(
      mean_val   = mean(value, na.rm = TRUE),
      max_val    = max(value, na.rm = TRUE),
      min_val    = min(value, na.rm = TRUE),
      cellsd_val = mean(cellsd, na.rm = TRUE),  # typical within-day cross-cell spread
      .groups = "drop"
    )

  daily_inst <- inst |>
    pivot_wider(id_cols = date, names_from = code,
                values_from = c(mean_val, max_val, min_val, cellsd_val),
                names_glue = "{code}_{.value}")

  next00 <- get_next_day_00h_reservoir(year, month)
  acc <- forcing |> filter(kind == "accumulated")
  days_in_month <- unique(as.Date(acc$time))

  acc_daily <- purrr::map_dfr(days_in_month, function(d) {
    next_d <- d + 1
    if (next_d %in% as.Date(acc$time)) {
      row <- acc |> filter(as.Date(time) == next_d, format(time, "%H:%M") == "00:00") |>
        select(code, value, cellsd)
    } else if (!is.null(next00)) {
      row <- next00 |> select(code, value = next00, cellsd = next00_cellsd)
    } else {
      return(tibble(date = d))   # end of record, no next-day value available
    }
    row |> mutate(date = d)
  })

  # BUG FOUND AND FIXED (3 Oct 2026, after the first full 1993-2025 run
  # produced two bogus 100%-NA columns "NA_dayTotal"/"NA_dayTotal_cellsd"):
  # the end-of-record branch above (line ~349, `tibble(date = d)`, no
  # `code` column) gets row-bound by map_dfr alongside rows that DO have
  # a `code` column - map_dfr fills the missing column with a literal
  # NA rather than leaving it absent, and that NA `code` then survives
  # into `pivot_wider()`, which dutifully creates a column named
  # "NA_dayTotal" from it. build_daily_forcing.R never hit this because
  # it pivots PER DAY inside the map function, before any end-of-record
  # row ever gets mixed in with a real one; this script pivots once for
  # the whole month instead, so the structurally-different row leaks
  # through. Fix: drop the code-NA row before pivoting - this loses
  # nothing, since every real variable already correctly comes out NA
  # for that date anyway (pivot_wider fills any date x code combination
  # that never appeared in the long data with NA - that's exactly the
  # documented, expected single end-of-record gap, untouched by this
  # fix). Confirmed against the already-written
  # daily_forcing_reservoir_1993-2025.csv: every real *_dayTotal column
  # was already correctly NA-at-2025-12-31-only; only the two synthetic
  # NA-named columns were ever wrong, and both were 100% NA throughout -
  # no real data was corrupted, just two dead columns to drop.
  acc_daily <- acc_daily |> filter(!is.na(code))

  acc_wide_val <- acc_daily |>
    pivot_wider(id_cols = date, names_from = code, values_from = value,
                names_glue = "{code}_dayTotal")
  acc_wide_sd <- if ("cellsd" %in% names(acc_daily)) {
    acc_daily |> select(date, code, cellsd) |>
      pivot_wider(id_cols = date, names_from = code, values_from = cellsd,
                  names_glue = "{code}_dayTotal_cellsd")
  } else {
    tibble(date = days_in_month)
  }

  full_join(daily_inst, acc_wide_val, by = "date") |> full_join(acc_wide_sd, by = "date")
}

# ---------------------------------------------------------------
# RUN_YEARS - 1993 test run reviewed and accepted 2 Oct 2026 (see
# Ataturk_Master_Checklist_v2.md, Stage C.2): W_JRC/W_CL diagnostic
# grids matched expectations exactly, 365 rows written, u2 moved into
# the expected 1.5-3 m/s range (1.671 -> 1.897, +0.226), licd dropped
# sharply (0.017 -> 0.005 mean, 0.144 -> 0.044 peak), and the t2m/ssrd
# sanity checks behaved as expected once the larger t2m shift was
# traced to a physically sound cause (see checklist). Widened to the
# full record.
# ---------------------------------------------------------------
RUN_YEARS <- 1993:2025

all_months <- expand.grid(month = 1:12, year = RUN_YEARS)
t_start <- Sys.time()
cat(sprintf("Starting %d months at %s...\n", nrow(all_months), format(t_start, "%H:%M:%S")))
result <- purrr::map2_dfr(all_months$year, all_months$month, build_month_reservoir) |>
  arrange(date)
t_end <- Sys.time()
cat(sprintf("Finished all months at %s (elapsed: %s)\n",
            format(t_end, "%H:%M:%S"), format(round(t_end - t_start, 2))))

out_path <- file.path(out_dir, sprintf("daily_forcing_reservoir_%s.csv",
                                        paste(range(RUN_YEARS), collapse = "-")))
write.csv(result, out_path, row.names = FALSE)
cat("Wrote", nrow(result), "reservoir-weighted daily rows to", out_path, "\n")
cat("\nColumn names:\n"); print(names(result))

# --- comparison vs the existing box-mean file - the whole point -----
cat("\n=== Comparison: reservoir-weighted vs box-mean, same date range ===\n")
box_mean_path <- here::here("data", "daily_forcing_1993-2025.csv")
if (file.exists(box_mean_path)) {
  bm <- read.csv(box_mean_path) |> mutate(date = as.Date(date)) |> filter(date %in% result$date)

  compare_var <- function(col, label) {
    if (!(col %in% names(bm)) || !(col %in% names(result))) {
      cat(sprintf("%-20s column '%s' missing from one file - skipped\n", label, col)); return(invisible())
    }
    merged <- inner_join(bm |> select(date, box_val = all_of(col)),
                          result |> select(date, res_val = all_of(col)), by = "date")
    cat(sprintf(
      "%-20s box-mean: %7.3f (%.3f to %.3f)   reservoir-mean: %7.3f (%.3f to %.3f)   mean diff: %+.3f\n",
      label,
      mean(merged$box_val, na.rm = TRUE), min(merged$box_val, na.rm = TRUE), max(merged$box_val, na.rm = TRUE),
      mean(merged$res_val, na.rm = TRUE), min(merged$res_val, na.rm = TRUE), max(merged$res_val, na.rm = TRUE),
      mean(merged$res_val - merged$box_val, na.rm = TRUE)))
  }

  compare_var("wspd10_mean_val", "u2 (wind, m/s)")
  compare_var("licd_mean_val",   "lake ice depth (m)")
  compare_var("t2m_mean_val",    "t2m (C) [sanity: should barely move]")
  compare_var("ssrd_dayTotal",   "ssrd (MJ/m2/day) [sanity: should barely move]")
} else {
  cat("Box-mean file not found at", box_mean_path, "- skipping comparison.\n")
}

cat("\n=== Between-cell spread (reservoir cells only) - overall summary ===\n")
spread_cols <- grep("_cellsd", names(result), value = TRUE)
for (col in spread_cols) {
  v <- result[[col]]
  cat(sprintf("%-28s mean=%.4f  range=%.4f to %.4f  (NA count: %d)\n",
              col, mean(v, na.rm = TRUE), min(v, na.rm = TRUE), max(v, na.rm = TRUE), sum(is.na(v))))
}
cat("\nLarge spread for a variable means the ~34 overlapping cells disagree a lot\n")
cat("with each other that day/on that layer - worth a look if any of these are\n")
cat("surprisingly large relative to the variable's own mean.\n")
