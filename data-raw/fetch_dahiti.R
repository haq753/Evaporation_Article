# ===============================================================
# data-raw/fetch_dahiti.R
#
# Stage C.2 (Ataturk_Master_Checklist_v2.md): DAHITI data for Atatürk
# Reservoir (DAHITI target ID 112, per the checklist - DAHITI's own ID
# system, unrelated to any ERA5-Land/CDS identifier).
#
# --- 1 Oct 2026: history, condensed after resolution -------------------
# This script went through three real iterations before landing here:
#   1. First version used DAHITI's v1 endpoint (username+password POST
#      fields) - got HTTP 404 (HTML homepage back) on every action.
#      DAHITI had moved to API v2 since that version was written:
#      single api_key auth (not username/password), one versioned
#      endpoint per action (`/api/v2/download-<action>/`) instead of
#      one endpoint + an `action` field, and different field names
#      (confirmed directly from DAHITI's own v2 doc pages, one page per
#      action, 1 Oct 2026).
#   2. Second version (v2 endpoints, api_key auth) ran against Aziz's
#      real account and got a genuinely mixed result:
#        - water_level: worked, but only after adding format=json
#          explicitly - the docs claim 'json' is the default, but the
#          real response without it was DAHITI's plain-text ascii
#          export, not JSON. (DAHITI's own example code always sets
#          format=json explicitly too - the default isn't trustworthy
#          in practice.) That ascii response's header also independently
#          confirmed DAHITI_ID 112 = Atatürk, Turkey, lon 38.5648 / lat
#          37.5794 - matching the Wikipedia-verified reservoir location.
#        - surface_area, volume_variation, hypsometry,
#          water_level_hypsometry: all HTTP 403 "Permission Denied".
#   3. **RESOLVED, not a script problem**: Aziz checked DAHITI's own
#      website for target 112 directly and confirmed via screenshot that
#      Atatürk's page shows only the "Water Level" product badge active
#      (blue) - Surface Area, Water Occurrence Mask, Land-Water Masks,
#      Volume Variations, Hypsometry, Bathymetry, and Water Level
#      (Hypsometry) are all greyed out, meaning DAHITI has never
#      computed those products for this reservoir, for anyone. The
#      "Permission Denied" wording was DAHITI's (confusing) way of
#      saying "no such product exists for this target" - not an
#      account/API-key restriction, and nothing to request from their
#      support. This script now only requests water_level - repeatedly
#      hitting four endpoints confirmed permanently unavailable for this
#      target would just waste time and clutter the output on every
#      future re-run.
#
# CONSEQUENCE FOR THE PAPER (see Ataturk_Master_Checklist_v2.md, Stage
# C.2): DAHITI contributes water-level altimetry only for this
# reservoir. It has no surface-area or volume product to reconcile
# against JRC - the "DAHITI vs. JRC area" and "847 vs. ~817 km^2"
# checklist items must be settled from JRC and the literature alone.
#
# CREDENTIALS: put your API key in .Renviron (Sys.getenv() reads it
# automatically; .Renviron is already outside git per this project's
# existing practice for secrets):
#   DAHITI_API_KEY=your_64_char_hex_key_here
# Generate it at https://dahiti.dgfi.tum.de -> My Profile -> API-Key.
# ===============================================================

library(httr)
library(jsonlite)
library(dplyr)

DAHITI_ID <- 112   # Ataturk Reservoir - independently confirmed above
                   # via DAHITI's own returned target metadata, not just
                   # assumed correct from the checklist.

api_key <- Sys.getenv("DAHITI_API_KEY")
if (api_key == "") {
  stop("DAHITI_API_KEY not set. Add it to .Renviron (see the comment block ",
       "at the top of this file) and restart R. Generate the key at ",
       "https://dahiti.dgfi.tum.de -> My Profile -> API-Key.")
}

BASE_URL <- "https://dahiti.dgfi.tum.de/api/v2/"
out_dir <- here::here("data")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

`%||%` <- function(a, b) if (is.null(a)) b else a

# get_field() exists because jsonlite doesn't always hand back a plain
# named list for a flat JSON object of all-scalar fields - it can come
# back as a named atomic vector instead (confirmed 1 Oct 2026: this
# crashed with "$ operator is invalid for atomic vectors" on the real
# response). This works for either shape, plus data.frame-shaped
# targets, instead of assuming one.
get_field <- function(x, field) {
  if (is.null(x)) return(NULL)
  if (is.list(x)) return(x[[field]])
  if (!is.null(names(x)) && field %in% names(x)) return(unname(x[[field]]))
  NULL
}

# --- fetch water_level - the only product DAHITI has for this target --
url <- paste0(BASE_URL, "download-water-level/")
resp <- POST(
  url,
  # format = "json" explicit, not relied on as a default - see the
  # dated note above: this came back as plain-text ascii export without it.
  body = list(api_key = api_key, dahiti_id = DAHITI_ID, format = "json"),
  encode = "form"
)

if (http_error(resp)) {
  stop("DAHITI request failed: HTTP ", status_code(resp), "\nRaw response:\n",
       content(resp, "text", encoding = "UTF-8"))
}

raw_text <- content(resp, "text", encoding = "UTF-8")
parsed <- tryCatch(
  fromJSON(raw_text),
  error = function(e) {
    stop("DAHITI response was not valid JSON. Raw response:\n", raw_text)
  }
)

if (!is.null(parsed$code) && parsed$code != 200) {
  stop("DAHITI reported an error: code ", parsed$code,
       ", message: ", parsed$message %||% "(none given)")
}

df <- as.data.frame(parsed$data)
out_path <- file.path(out_dir, "dahiti_water_level.csv")
write.csv(df, out_path, row.names = FALSE)
cat("OK - ", nrow(df), " rows -> ", out_path, "\n", sep = "")
cat("Columns:", paste(names(df), collapse = ", "), "\n")
# The documented schema said "date"; the real response's data column
# turned out to be "datetime" instead (confirmed 1 Oct 2026 - another
# case of the docs not matching the live API exactly) - check both
# rather than hardcoding one and silently reporting nothing.
date_col <- intersect(c("date", "datetime"), names(df))[1]
if (!is.na(date_col)) {
  cat("Coverage:", min(df[[date_col]]), "to", max(df[[date_col]]), "\n")
  cat("(Checklist item 'note coverage gaps, especially pre-2010' - check the\n")
  cat("actual date spacing above, don't assume even coverage across the range.)\n")
} else {
  cat("(no date/datetime column found - see Columns above)\n")
}

# --- sanity check: does DAHITI ID 112 actually mean Ataturk? --------
# The field names assumed here (target_name/country/longitude/latitude)
# were taken from DAHITI's documented example schema, which just proved
# unreliable for the data columns too (date vs datetime) - rather than
# guess at another set of assumed names and risk the same silent
# "(none)" failure again, this prints the RAW structure so the actual
# field names are visible directly, however DAHITI is really shaping
# this response right now.
cat("\n=== Target ID sanity check (is DAHITI ID", DAHITI_ID, "really Ataturk?) ===\n")
target <- parsed$target
if (!is.null(target)) {
  cat("Raw target object as returned by DAHITI (field names may differ from\n")
  cat("the docs - that's exactly what broke last time, so showing the real\n")
  cat("thing instead of guessing at names again):\n")
  print(target)
  cat("\n(Compare the values above by eye against the Wikipedia-verified Ataturk\n")
  cat("location: Turkey, lon ~38.3-38.7E, lat ~37.4-37.6N - 539 points/rows\n")
  cat("already matches DAHITI's own prior count for this target, which is\n")
  cat("already strong evidence ID 112 is right regardless of this printout.)\n")
} else {
  cat("  No 'target' key in this response at all - top-level keys were:",
      paste(names(parsed), collapse = ", "), "\n")
}

cat("\nNOTE: surface_area, volume_variation, hypsometry, and water_level_hypsometry\n")
cat("are not requested here - confirmed 1 Oct 2026 (DAHITI's own website for target\n")
cat("112) that none of those products exist for Ataturk. Water level is the only\n")
cat("DAHITI product available for this reservoir; see the checklist for what that\n")
cat("means for the 'DAHITI vs JRC area' and '847 vs 817 km^2' items.\n")
