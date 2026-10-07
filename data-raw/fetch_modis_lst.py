"""
data-raw/fetch_modis_lst.py

Stage C.3 (Ataturk_Master_Checklist_v2.md): MODIS Land Surface Temperature
for Atatürk Reservoir, via Earth Engine. Produces:

  1. `data/modis_lst_2000_2025.csv` - one row per 8-day composite per
     product (Terra/MOD11A2, Aqua/MYD11A2), with QA-masked reservoir-
     mean day and night LST (Celsius) plus several QA diagnostic
     fractions (see below). Kept as two separate rows per date (one per
     satellite) rather than pre-merged into one series - which satellite
     to trust/average on any given date is a decision for later in the
     paper, not baked in here.
  2. `data/modis_lst_cloud_contamination_by_month.csv` - the checklist's
     explicit "report cloud-contamination rate by month" item, derived
     from the same per-image data (mean cloud fraction per calendar
     month, across all years, day and night, both satellites).

Dataset facts below were verified against Google's live Earth Engine
catalog pages and the authoritative LP DAAC MOD11 Collection 6.1 User
Guide (3-4 Oct 2026 research), not assumed or recalled from training -
same discipline as fetch_jrc_surface_water.py:

  - Collection IDs confirmed from the catalog, not guessed from the
    checklist's shorthand "MOD11A2 + MYD11A2":
      Terra:  MODIS/061/MOD11A2  (v6.1 - the current, non-deprecated
              version; there is no newer collection as of this check)
      Aqua:   MODIS/061/MYD11A2
    (v6.1 replaced v6.0, which Google's catalog marks superseded - the
    same "check which version is current, don't assume the checklist's
    implicit version is still right" step as the JRC v1.4-vs-v1.3 case.)
  - Temporal coverage, confirmed from the catalog: MOD11A2.061 starts
    2000-02-18, MYD11A2.061 starts 2002-07-04; both show as actively
    updated through at least Sept 2026 - i.e. BOTH products genuinely
    cover the full requested 2000-2025 window (Aqua necessarily misses
    2000-02-18 to 2002-07-03, which is expected and fine: Terra alone
    covers that gap).
  - LST_Day_1km / LST_Night_1km: scale factor 0.02, units Kelvin -
    checklist's "apply scale factor 0.02" is confirmed CORRECT, not
    stale. Converted to Celsius here (-273.15) to match every other
    temperature column already in this project's forcing data.
  - QC_Day / QC_Night bit layout, confirmed against the LP DAAC MOD11
    Collection 6.1 User Guide (Table 18), not just the checklist's
    shorthand:
      bits 0-1 (mandatory QA):     00 = LST produced, good quality
                                    01 = LST produced, other quality
                                         (examine more detailed QA)
                                    10 = LST NOT produced, cloud effects
                                    11 = LST NOT produced, other reason
      bits 2-3 (data quality):     00 = good, 01 = other, 10/11 = TBD
      bits 4-5 (emissivity error): 00 <=0.01, 01 <=0.02, 10 <=0.04,
                                    11 > 0.04
      bits 6-7 (LST error):        00 <=1K, 01 <=2K, 10 <=3K, 11 > 3K

  *** QA MASK DESIGN - CHANGED after the first run (3-4 Oct 2026) ***
  The first version of this script kept ONLY pixels with mandatory QA
  == 00 ("good quality") and reported that as the LST mean. Running it
  over the reservoir returned a *0% acceptance rate across all 2,269
  composites, 2000-2025* - not a handful of cloudy months, literally
  every single 8-day window - while the independently-tracked cloud
  fraction (mandatory QA == 10) stayed low and seasonally sensible
  (0.3-5%). That pattern (near-zero "good quality" over water, cloud
  fraction normal) is a documented, known MODIS LST behaviour, not a
  bug in this script or the reservoir geometry:
    - The LP DAAC user guide itself gives water bodies a MORE relaxed
      cloud-confidence threshold for retrieval (clear-sky confidence
      >=66% over lakes vs. a stricter threshold over land), implying
      the algorithm already expects water pixels to behave differently
      from land pixels under the same QA scheme.
    - A published lake-LST validation study ("A Strict Validation of
      MODIS Lake Surface Water Temperature on the Tibetan Plateau",
      Remote Sensing / MDPI 2022) explicitly filters MODIS LST over
      lakes using mandatory-QA "good quality" **combined with an
      LST-error threshold of <=2K**, rather than mandatory-QA alone -
      i.e. published lake-LST work does not treat "good quality"
      (bits0-1==00) as the only usable class; it treats "LST produced"
      (bits0-1 in {00,01}) plus an explicit error-bit check as the
      usable set.
  Rather than guess a threshold, this script follows that published
  approach: accept a pixel when (a) LST was actually produced - i.e.
  mandatory QA is 00 or 01, NOT the two "not produced" classes - AND
  (b) the LST-error bits indicate <=2K average error (bits6-7 in
  {00,01}). This is a methods decision, not a re-derivable fact, so the
  2K threshold and the "any produced-with-low-error pixel is usable"
  choice should be stated explicitly in the paper's methods section
  when this data is used, and are flagged here for Aziz/Tombul to
  confirm rather than silently baked in. To make that check possible,
  the output keeps BOTH the new accepted-pixel fraction AND the original
  strict "good quality" (bits0-1==00 exactly) fraction as separate
  columns, so the 0%-vs-nonzero contrast that motivated this change is
  visible in the data itself, not just this comment.

  *** QUERY CHUNKING - added after the second run (4 Oct 2026) ***
  The Terra fetch (1,189 composites, one single .getInfo() call over
  the whole 2000-2025 collection) succeeded; the very next Aqua fetch
  (1,080 composites - FEWER images, same per-image reduceRegion cost)
  failed with "Computation timed out." Checked against Google's own
  Earth Engine documentation rather than just adding a blind retry:
  the "interactive" environment that a plain .getInfo() call runs in
  has a hard, documented 5-minute ceiling per request (Processing
  Environments guide), and Earth Engine's own debugging guide's fix for
  this specific error is to move long computations out of that
  synchronous path. A full asynchronous Export.table.toDrive() pipeline
  (task polling, Drive download) would work but is a much bigger change
  for this project than the problem needs. Instead, each product's
  collection is now fetched one calendar year at a time - 26 small
  .getInfo() calls instead of 1 large one - which keeps every single
  request comfortably under the 5-minute ceiling regardless of server-
  side load, and means a transient failure only costs one year's retry,
  not the whole collection. A short retry-with-pause is also added
  specifically for "timed out" errors, since a community report of the
  same error (geemap discussions) describes it recurring intermittently
  under load even after splitting into smaller requests - i.e. this is
  a known-flaky failure mode, not just a size problem, so a bare retry
  is worth having in addition to chunking.

PREREQUISITES: same Earth Engine Python environment already set up and
verified for fetch_jrc_surface_water.py (pip install earthengine-api;
ee.Authenticate() once; EE_PROJECT below reuses that same confirmed-
working Google Cloud project), PLUS shapely, which fetch_jrc_surface_
water.py didn't need (pip install shapely) - used here to re-apply the
same make_valid() fix to the exported polygon that the W_JRC weight
computation in remask_forcing_to_reservoir.R's header documents, since
the polygon file itself still has that one confirmed-negligible ring
self-intersection and Earth Engine's geometry constructor can be
strict about accepting invalid rings.

NOT done here (checklist marks it optional): Landsat thermal back to
1993 for the pre-MODIS period. Left for later - MODIS alone already
covers 2000-2025, which is the bulk of the 33-year record, and Landsat
thermal needs its own emissivity/atmospheric-correction verification
pass before it's worth building.
"""

import json
import os
import time
import ee
from shapely.geometry import shape, mapping
from shapely.validation import make_valid

# ---------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------
EE_PROJECT = "ataturk-evaporation"  # same confirmed-working project as fetch_jrc_surface_water.py

DATE_START = "2000-01-01"
DATE_END = "2026-01-01"  # exclusive upper bound - covers all of 2025

# LST-error bit class accepted alongside "LST produced" mandatory QA
# (bits0-1 in {0,1}) - see the QA MASK DESIGN note above. Class 1
# means "average LST error <= 2K" (bits6-7 == 00 or 01). Chosen to
# match the published lake-LST precedent (Tibetan Plateau lake study);
# revisit if Aziz/Tombul want a tighter (0 -> <=1K) or looser threshold
# for the paper's stated methodology.
ACCEPT_MAX_LST_ERROR_CLASS = 1

# Retry behaviour for Earth Engine's "Computation timed out" error - see
# the QUERY CHUNKING note above. A retry only fires for that specific
# error, never for a genuine data/auth problem.
GETINFO_MAX_ATTEMPTS = 3
GETINFO_RETRY_WAIT_S = 15

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "..", "data")
os.makedirs(DATA_DIR, exist_ok=True)

POLYGON_GEOJSON_PATH = os.path.join(DATA_DIR, "ataturk_jrc_max_extent_polygon.geojson")
OUT_LST_CSV = os.path.join(DATA_DIR, "modis_lst_2000_2025.csv")
OUT_CLOUD_BY_MONTH_CSV = os.path.join(DATA_DIR, "modis_lst_cloud_contamination_by_month.csv")

# Same reservoir definition as the ERA5-Land re-mask (Stage C.2,
# remask_forcing_to_reservoir.R) - confirmed with Aziz 1 Oct 2026,
# original box. Using the SAME polygon file here, not re-deriving it,
# keeps "the reservoir" consistent across every extraction in the paper.
LST_SCALE_M = 1000  # native MODIS LST resolution


def load_reservoir_geometry():
    """Load the final reservoir polygon and fix it the same way the
    W_JRC weight computation did (remask_forcing_to_reservoir.R's
    header documents this exact fix) - the exported GeoJSON has one
    confirmed-negligible ring self-intersection (~9e-16 deg^2), which
    Earth Engine's geometry constructor can be strict about. Fixed at
    load time here rather than hand-editing the file, so this script
    always tracks whatever the authoritative polygon file currently
    contains."""
    with open(POLYGON_GEOJSON_PATH, encoding="utf-8") as f:
        gj = json.load(f)
    raw_geom = gj["geometry"] if gj.get("type") == "Feature" else gj
    fixed = make_valid(shape(raw_geom))
    fixed_geojson = mapping(fixed)
    print(f"Reservoir polygon loaded and validated: geometry type = {fixed_geojson['type']}")
    return ee.Geometry(fixed_geojson)


def build_reducer_fn(reservoir_geom):
    """Per-image reducer. For each of day/night:
      - `strict_good_*_frac`: fraction of pixels with mandatory QA
        bits0-1 == 00 exactly ("good quality") - kept as a diagnostic
        only, not used for the LST mean (see QA MASK DESIGN note).
      - `cloud_*_frac`: fraction with mandatory QA bits0-1 == 10
        ("not produced, cloud") - the checklist's cloud-contamination
        metric, unchanged by the QA-mask redesign.
      - `accept_*_frac`: fraction actually used for the LST mean -
        mandatory QA bits0-1 in {00,01} ("LST produced", confidence
        either way) AND LST-error bits6-7 <= ACCEPT_MAX_LST_ERROR_CLASS.
      - `lst_*_c`: LST mean (Celsius) over only the accepted pixels.
    Combined into one reduceRegion call per image rather than many
    separate calls, to keep the per-image server-side cost down."""

    def reduce_image(img):
        qc_day = img.select("QC_Day")
        qc_night = img.select("QC_Night")

        mandatory_day = qc_day.bitwiseAnd(3)              # bits 0-1
        mandatory_night = qc_night.bitwiseAnd(3)
        lst_err_day = qc_day.rightShift(6).bitwiseAnd(3)    # bits 6-7
        lst_err_night = qc_night.rightShift(6).bitwiseAnd(3)

        strict_good_day = mandatory_day.eq(0)
        strict_good_night = mandatory_night.eq(0)
        cloud_day = mandatory_day.eq(2)
        cloud_night = mandatory_night.eq(2)

        produced_day = mandatory_day.lt(2)    # 00 or 01 -> "LST produced"
        produced_night = mandatory_night.lt(2)
        accept_day = produced_day.And(lst_err_day.lte(ACCEPT_MAX_LST_ERROR_CLASS))
        accept_night = produced_night.And(lst_err_night.lte(ACCEPT_MAX_LST_ERROR_CLASS))

        lst_day_c = (
            img.select("LST_Day_1km").multiply(0.02).subtract(273.15)
            .updateMask(accept_day).rename("lst_day_c")
        )
        lst_night_c = (
            img.select("LST_Night_1km").multiply(0.02).subtract(273.15)
            .updateMask(accept_night).rename("lst_night_c")
        )

        stack = ee.Image.cat([
            lst_day_c,
            lst_night_c,
            accept_day.rename("accept_day_frac"),
            accept_night.rename("accept_night_frac"),
            strict_good_day.rename("strict_good_day_frac"),
            strict_good_night.rename("strict_good_night_frac"),
            cloud_day.rename("cloud_day_frac"),
            cloud_night.rename("cloud_night_frac"),
        ])

        stats = stack.reduceRegion(
            reducer=ee.Reducer.mean(),
            geometry=reservoir_geom,
            scale=LST_SCALE_M,
            maxPixels=1e9,
        )

        return ee.Feature(None, stats).set(
            "date", img.date().format("YYYY-MM-dd")
        )

    return reduce_image


def _getinfo_with_retry(fc, label):
    """getInfo() on a mapped FeatureCollection, retrying only on Earth
    Engine's "Computation timed out" error - see the QUERY CHUNKING
    note in the module docstring. Any other exception (auth, quota,
    genuine bug) is raised immediately, not retried."""
    last_err = None
    for attempt in range(1, GETINFO_MAX_ATTEMPTS + 1):
        try:
            return fc.getInfo()
        except ee.ee_exception.EEException as e:
            last_err = e
            if "timed out" in str(e).lower() and attempt < GETINFO_MAX_ATTEMPTS:
                print(f"    {label}: Earth Engine computation timed out "
                      f"(attempt {attempt}/{GETINFO_MAX_ATTEMPTS}) - this is the documented "
                      f"5-minute interactive-request ceiling, not specific to this chunk; "
                      f"waiting {GETINFO_RETRY_WAIT_S}s and retrying the same request.")
                time.sleep(GETINFO_RETRY_WAIT_S)
                continue
            raise
    raise last_err  # pragma: no cover - loop above always returns or raises


def fetch_product(collection_id, label, reservoir_geom):
    """Fetch one product's whole date range, one calendar year at a
    time (see QUERY CHUNKING note above) - 26 small .getInfo() calls
    instead of one large one that can exceed Earth Engine's 5-minute
    interactive-request limit."""
    base_coll = ee.ImageCollection(collection_id).filterBounds(reservoir_geom)

    start_year = int(DATE_START[:4])
    end_year = int(DATE_END[:4])  # DATE_END is exclusive, e.g. "2026-01-01" -> loop ends at 2025

    rows = []
    total_n = 0
    print(f"\n{label} ({collection_id}): fetching {start_year}-{end_year - 1} year by year")
    for year in range(start_year, end_year):
        chunk_start = max(f"{year}-01-01", DATE_START)
        chunk_end = min(f"{year + 1}-01-01", DATE_END)
        if chunk_start >= chunk_end:
            continue

        coll_year = base_coll.filterDate(chunk_start, chunk_end)
        n = coll_year.size().getInfo()
        if n == 0:
            print(f"  {year}: 0 composites, skipped")
            continue

        fc = coll_year.map(build_reducer_fn(reservoir_geom))
        result = _getinfo_with_retry(fc, label=f"{label} {year}")

        for feat in result["features"]:
            p = feat["properties"]
            rows.append({
                "date": p.get("date"),
                "product": label,
                "lst_day_c": p.get("lst_day_c"),
                "lst_night_c": p.get("lst_night_c"),
                "accept_day_frac": p.get("accept_day_frac"),
                "accept_night_frac": p.get("accept_night_frac"),
                "strict_good_day_frac": p.get("strict_good_day_frac"),
                "strict_good_night_frac": p.get("strict_good_night_frac"),
                "cloud_day_frac": p.get("cloud_day_frac"),
                "cloud_night_frac": p.get("cloud_night_frac"),
            })
        total_n += n
        print(f"  {year}: {n} composites (cumulative {total_n})")

    print(f"{label}: {total_n} total 8-day composites fetched")
    return rows


def main():
    t_start = time.time()
    ee.Authenticate()  # no-op if already authenticated on this machine
    ee.Initialize(project=EE_PROJECT)

    reservoir_geom = load_reservoir_geometry()

    rows = []
    rows += fetch_product("MODIS/061/MOD11A2", "terra", reservoir_geom)
    rows += fetch_product("MODIS/061/MYD11A2", "aqua", reservoir_geom)
    rows.sort(key=lambda r: (r["date"], r["product"]))

    # --- write the main LST CSV --------------------------------
    fieldnames = ["date", "product", "lst_day_c", "lst_night_c",
                  "accept_day_frac", "accept_night_frac",
                  "strict_good_day_frac", "strict_good_night_frac",
                  "cloud_day_frac", "cloud_night_frac"]
    import csv
    with open(OUT_LST_CSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"\nWrote {len(rows)} rows to {OUT_LST_CSV}")

    n_missing_lst = sum(1 for r in rows if r["lst_day_c"] is None or r["lst_night_c"] is None)
    print(f"Rows with no accepted (QC-produced, LST-error<=2K) pixel for day OR night: "
          f"{n_missing_lst} ({100*n_missing_lst/len(rows):.1f}%) - expected to be occasionally "
          "nonzero (fully cloud/no-retrieval 8-day windows happen), not a bug by itself; look "
          "at the by-month breakdown below for whether it clusters seasonally.")

    # --- whole-record QA diagnostic: accepted vs strict-good vs cloud --
    def _mean(key):
        vals = [r[key] for r in rows if r[key] is not None]
        return sum(vals) / len(vals) if vals else float("nan")

    print("\n=== QA diagnostic, whole record (confirms/contradicts the QA MASK DESIGN note) ===")
    print(f"  day:   accepted={_mean('accept_day_frac'):.3f}  "
          f"strict-good(QC==0 exactly)={_mean('strict_good_day_frac'):.3f}  "
          f"cloud={_mean('cloud_day_frac'):.3f}")
    print(f"  night: accepted={_mean('accept_night_frac'):.3f}  "
          f"strict-good(QC==0 exactly)={_mean('strict_good_night_frac'):.3f}  "
          f"cloud={_mean('cloud_night_frac'):.3f}")

    # --- cloud-contamination rate by calendar month --------------
    from collections import defaultdict
    by_month = defaultdict(lambda: {"cloud_day": [], "cloud_night": []})
    for r in rows:
        month = int(r["date"][5:7])
        if r["cloud_day_frac"] is not None:
            by_month[month]["cloud_day"].append(r["cloud_day_frac"])
        if r["cloud_night_frac"] is not None:
            by_month[month]["cloud_night"].append(r["cloud_night_frac"])

    with open(OUT_CLOUD_BY_MONTH_CSV, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["month", "mean_cloud_day_frac", "mean_cloud_night_frac", "n_composites"])
        print("\n=== Cloud-contamination rate by calendar month (both satellites combined) ===")
        for month in range(1, 13):
            day_vals = by_month[month]["cloud_day"]
            night_vals = by_month[month]["cloud_night"]
            mean_day = sum(day_vals) / len(day_vals) if day_vals else float("nan")
            mean_night = sum(night_vals) / len(night_vals) if night_vals else float("nan")
            writer.writerow([month, mean_day, mean_night, len(day_vals)])
            print(f"  month {month:2d}: day={mean_day:.3f}  night={mean_night:.3f}  "
                  f"(n={len(day_vals)} composites)")
    print(f"\nWrote cloud-contamination-by-month table to {OUT_CLOUD_BY_MONTH_CSV}")
    print("\nIf any month's day/night fraction is dramatically higher than the rest, "
          "that's the seasonal cloud pattern this item asked to check for - plausible "
          "candidates would be spring/early-summer convective cloud in this region, but "
          "look at the actual numbers rather than assume that's what shows up.")

    t_end = time.time()
    print(f"\nTotal run time: {(t_end - t_start) / 60:.1f} minutes")


if __name__ == "__main__":
    main()
