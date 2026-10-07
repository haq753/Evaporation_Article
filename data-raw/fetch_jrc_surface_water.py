"""
data-raw/fetch_jrc_surface_water.py

Stage C.2 (Ataturk_Master_Checklist_v2.md): JRC Global Surface Water via
Earth Engine. Two outputs, both needed downstream:

  1. A reservoir polygon (GeoJSON) built from the `max_extent` band of
     JRC/GSW1_4/GlobalSurfaceWater - "this becomes the mask for every
     other extraction" per the checklist, so it needs to be right before
     Stage D touches it.
  2. A monthly reservoir surface-area time series from
     JRC/GSW1_4/MonthlyHistory, for direct comparison against DAHITI's
     own surface-area product (the "reconcile DAHITI vs. JRC area" and
     "resolve the 847 vs. ~817 km^2 question" checklist items).

Dataset facts below were verified against the live Earth Engine catalog
pages before writing this (29 Sep-1 Oct 2026 research), not assumed:

  - JRC/GSW1_4/GlobalSurfaceWater (v1.4) is the CURRENT version - v1.3
    and v1.0 are explicitly marked "deprecated" on Google's own catalog,
    v1.4 is not, and there is no v1.5 yet. The checklist's dataset ID is
    correct and current.
  - GlobalSurfaceWater bands used here:
      max_extent  - binary, 1 = water detected at ANY point 1984-2021.
                    This is the basis for the reservoir polygon.
      occurrence  - % of time water was present (0-100). Used only as a
                    sensitivity check alongside max_extent, for the
                    847-vs-817 question - NOT as the polygon itself.
  - JRC/GSW1_4/MonthlyHistory: 454 monthly images, March 1984 - December
    2021 (confirmed from the catalog, not assumed). THIS MEANS THERE IS
    NO JRC COVERAGE FOR 2022-2025 - a real ~4-year gap against this
    study's 1993-2025 period, not something build_month-style code can
    paper over. DAHITI's altimetry-based product (pulled separately,
    see fetch_dahiti.R) needs to be checked for whether it extends past
    2021; if not, 2022-2025 area will need a documented
    extrapolation/assumption at Stage E, not silent interpolation.

--- 1 Oct 2026: box-edge investigation, RESOLVED - box reverted ---------
The first run of this script (original box, same as the completed
ERA5-Land download) reported "all four edges touch the box", via a naive
check on the raw max_extent image. That was a false positive on three
sides: by-hand vertex analysis of the exported polygon showed 0 of
44,625 vertices within 1 pixel of the west/south/north edges, and only
18 right on the east edge (39.05003 vs. box edge 39.05000, ~3 m).

Widening EAST to 39.5 and re-running pushed the area to 840.49 km^2
(now ~3% ABOVE the 817 km^2 literature figure) and made the NORTH edge
start touching too. Rather than keep widening (which never converges if
the thing beyond the edge is a river, not a lake), the exported polygon
was plotted directly (matplotlib, against the dam wall and reservoir
centroid from Wikipedia). The picture settled it: the polygon is a
broad, dendritic, many-armed lake out to roughly 39.0-39.1 E, then
becomes a long, narrow, uniformly-wide, winding channel continuing
northeast - visually and geomorphologically the Euphrates/tributary
channel upstream of the reservoir's real backwater limit, not reservoir
surface. That river is what was touching the east edge, and then the
north edge, as the box widened - confirmed NOT a missing reservoir lobe.

DECISION (confirmed with Aziz, 1 Oct 2026): keep the ORIGINAL box
`[38.1, 38.05, 37.35, 39.05]` - the same one already used for the
complete, 792-file ERA5-Land download - rather than widen further. The
~16 km^2 (~2%) gap between max_extent-in-this-box (801.36 km^2) and the
817 km^2 literature figure is accepted and documented as a minor,
transitional reservoir-to-river-arm effect, not chased further. No
ERA5-Land re-download is needed. This script is reverted to that
original box; the widened-box run's outputs (polygon + monthly CSV) are
superseded by this run's and should not be used.

Because the box is unchanged from the very first run, the east edge is
EXPECTED to still show as "touching" below - that's the accepted,
understood river-channel effect above, not a new problem. West/south/
north touching, if it ever happened, WOULD be new and worth
investigating, since those were verified clean against this exact box.

PREREQUISITES (one-time, on your machine, not something this script can
do for you):
  pip install earthengine-api
  python -c "import ee; ee.Authenticate()"   # opens a browser login once
  You also need a Google Cloud project with the Earth Engine API enabled
  - if you already use Earth Engine you likely have one; if not, create
  one free at https://console.cloud.google.com/ and enable the Earth
  Engine API for it, then put its project ID in EE_PROJECT below.
"""

import json
import os
import ee

# ---------------------------------------------------------------
# CONFIG
# ---------------------------------------------------------------
EE_PROJECT = "ataturk-evaporation"  # confirmed working project ID, 1 Oct 2026

# ORIGINAL box, same as the completed ERA5-Land download
# (data-raw/download_era5land_bulk.R) - reverted here after the box-edge
# investigation above confirmed widening further just chases the
# Euphrates channel upstream rather than finding more reservoir.
NORTH, WEST, SOUTH, EAST = 38.1, 38.05, 37.35, 39.05

# Paths relative to THIS FILE's location, not the working directory -
# R's here::here() does this automatically; Python doesn't, so it's
# done explicitly here. Assumes this file lives in Code/data-raw/,
# same as every other data-raw script in this project.
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "..", "data")
os.makedirs(DATA_DIR, exist_ok=True)

OUT_POLYGON_GEOJSON = os.path.join(DATA_DIR, "ataturk_jrc_max_extent_polygon.geojson")
OUT_MONTHLY_AREA_CSV = os.path.join(DATA_DIR, "jrc_monthly_surface_area_1984_2021.csv")

# ~1 pixel at 30 m resolution, in degrees - used to decide "is this
# vertex essentially AT the box edge" rather than requiring an exact
# floating-point match.
EDGE_TOL_DEG = 0.0003

# Sides that are EXPECTED to show as touching with the original box,
# per the investigation above (the Euphrates channel continuing past
# the reservoir's real limit) - printed as accepted, not alarming.
EXPECTED_TOUCH_SIDES = {"east"}


def main():
    ee.Authenticate()  # no-ops if already authenticated on this machine
    ee.Initialize(project=EE_PROJECT)

    aoi = ee.Geometry.Rectangle([WEST, SOUTH, EAST, NORTH])

    gsw = ee.Image("JRC/GSW1_4/GlobalSurfaceWater").clip(aoi)
    max_extent = gsw.select("max_extent")
    occurrence = gsw.select("occurrence")

    # --- sensitivity check for the 847 vs ~817 km^2 question ----
    # Four different definitions of "the reservoir area", computed from
    # the same source product, printed side by side rather than picking
    # one and hoping it matches. Whichever the DAHITI figure (pulled
    # separately) lines up with tells us which definition the
    # literature's 847/817 figures were probably using.
    def area_km2(binary_image, label):
        stat = binary_image.multiply(ee.Image.pixelArea()).reduceRegion(
            reducer=ee.Reducer.sum(), geometry=aoi, scale=30, maxPixels=1e10
        )
        # reduceRegion with pixelArea() already returns m^2 - convert here
        km2 = ee.Number(stat.values().get(0)).divide(1e6)
        print(f"{label}: {km2.getInfo():.2f} km^2")
        return km2

    print("=== Area sensitivity check (same source, different thresholds) ===")
    area_km2(max_extent, "max_extent (ever wet, 1984-2021)")
    area_km2(occurrence.gte(50), "occurrence >= 50%")
    area_km2(occurrence.gte(75), "occurrence >= 75%")
    area_km2(occurrence.gte(90), "occurrence >= 90% (near-permanent water)")

    # --- build the reservoir polygon ------------------------------
    vectors = max_extent.selfMask().reduceToVectors(
        geometry=aoi, scale=30, geometryType="polygon", maxPixels=1e10, eightConnected=True
    )
    # Keep only the largest polygon - drops any small unrelated ponds,
    # river threads, or irrigation canals elsewhere in the box that
    # max_extent also picked up (these are exactly what made the
    # earlier whole-image edge check misleading).
    vectors = vectors.map(lambda f: f.set("area_m2", f.geometry().area(1)))
    largest = vectors.sort("area_m2", False).first()
    largest_geom = largest.geometry()
    geojson = largest_geom.getInfo()

    with open(OUT_POLYGON_GEOJSON, "w") as f:
        json.dump({"type": "Feature", "geometry": geojson, "properties": {}}, f)
    print(f"\nWrote reservoir polygon to {OUT_POLYGON_GEOJSON}")
    print("This supersedes the widened-box run's polygon - that one included a stretch")
    print("of the Euphrates channel upstream of the reservoir (see the dated note above).")

    # --- edge-of-box check, on the ACTUAL RESERVOIR POLYGON ---------
    # Reads the exported polygon's own bounding coordinates and compares
    # them to each box edge directly, rather than testing the raw image
    # (which can't distinguish the reservoir from unrelated water
    # features elsewhere in the box - see the dated note above for why
    # that distinction mattered here).
    bounds_coords = largest_geom.bounds(1).coordinates().get(0).getInfo()
    lons = [c[0] for c in bounds_coords]
    lats = [c[1] for c in bounds_coords]
    poly_west, poly_east = min(lons), max(lons)
    poly_south, poly_north = min(lats), max(lats)

    print("\n=== Edge-of-box check (on the actual reservoir polygon) ===")
    print(f"Reservoir polygon bbox: west={poly_west:.5f} east={poly_east:.5f} "
          f"south={poly_south:.5f} north={poly_north:.5f}")
    print(f"AOI box:                west={WEST:.5f} east={EAST:.5f} "
          f"south={SOUTH:.5f} north={NORTH:.5f}")

    edges = {
        "west": (WEST - poly_west) if poly_west < WEST + EDGE_TOL_DEG else None,
        "east": (poly_east - EAST) if poly_east > EAST - EDGE_TOL_DEG else None,
        "south": (SOUTH - poly_south) if poly_south < SOUTH + EDGE_TOL_DEG else None,
        "north": (poly_north - NORTH) if poly_north > NORTH - EDGE_TOL_DEG else None,
    }
    unexpected_touch = False
    for side, margin_over in edges.items():
        if margin_over is not None:
            if side in EXPECTED_TOUCH_SIDES:
                print(f"  {side} edge: touches, as expected (Euphrates channel upstream of "
                      f"the reservoir, confirmed 1 Oct 2026 - not a missing-area problem)")
            else:
                unexpected_touch = True
                print(f"  {side} edge: TOUCHES and this is UNEXPECTED (this side was verified "
                      f"clean before) -> investigate before trusting the polygon")
        else:
            if side in ("west", "south"):
                margin = (poly_west if side == "west" else poly_south) - (WEST if side == "west" else SOUTH)
            else:
                margin = (EAST if side == "east" else NORTH) - (poly_east if side == "east" else poly_north)
            print(f"  {side} edge: clear, margin = {margin:.5f} deg (~{margin * 100:.1f} km)")

    if unexpected_touch:
        print("\nUNEXPECTED edge touch detected - this is new since the 1 Oct 2026 "
              "investigation. Do not trust this polygon as the mask until re-checked.")
    else:
        print("\nNo unexpected edge touches - matches the 1 Oct 2026 investigation "
              "(east touch is the known, accepted Euphrates-channel effect).")

    # --- monthly surface area time series --------------------------
    # NOTE: this computes a reduceRegion over 454 monthly images
    # synchronously via getInfo() - fine for a one-off, but if it times
    # out on your connection, the fallback is an async export instead:
    # ee.batch.Export.table.toDrive(collection=monthly.map(monthly_area),
    # description=..., fileFormat='CSV') and check the Earth Engine Code
    # Editor's Tasks tab for completion, rather than waiting on this
    # call - not implemented here since the synchronous path is simpler
    # and should work for a single reservoir-sized AOI.
    monthly = ee.ImageCollection("JRC/GSW1_4/MonthlyHistory").filterBounds(aoi)

    def monthly_area(img):
        water = img.eq(2).selfMask().clip(aoi)
        stat = water.multiply(ee.Image.pixelArea()).reduceRegion(
            reducer=ee.Reducer.sum(), geometry=aoi, scale=30, maxPixels=1e10
        )
        # .get(key, default): EE's reduceRegion can omit the "water" key
        # entirely when every pixel in the AOI is masked that month (e.g.
        # total cloud cover) - .get() WITHOUT a default raises a server-
        # side error in that case, which would silently kill the whole
        # collection.map() call rather than just flagging one bad month.
        # -1 m^2 is an unambiguous sentinel (a real area can't be
        # negative) carried through as -1 km^2 below, so "no data this
        # month" is visible in the output instead of looking like a
        # genuine 0 km^2 (dry reservoir) reading.
        water_m2 = ee.Number(stat.get("water", -1))
        area_km2_val = ee.Algorithms.If(water_m2.eq(-1), -1, water_m2.divide(1e6))
        return ee.Feature(
            None,
            {
                "year": img.get("year"),
                "month": img.get("month"),
                "area_km2": area_km2_val,
            },
        )

    records = monthly.map(monthly_area).getInfo()["features"]
    rows = [r["properties"] for r in records]
    rows.sort(key=lambda r: (r["year"], r["month"]))

    with open(OUT_MONTHLY_AREA_CSV, "w") as f:
        f.write("year,month,area_km2\n")
        for r in rows:
            f.write(f"{r['year']},{r['month']},{r['area_km2']}\n")

    n_no_data = sum(1 for r in rows if r["area_km2"] == -1)
    print(f"\nWrote {len(rows)} monthly rows to {OUT_MONTHLY_AREA_CSV}")
    print(f"  {n_no_data} month(s) had no valid pixels in the AOI (area_km2 = -1, not a real 0)")
    print(f"  Coverage: {rows[0]['year']}-{rows[0]['month']:02d} to {rows[-1]['year']}-{rows[-1]['month']:02d}")
    print("  REMINDER: this ends in Dec 2021 - there is no JRC coverage for 2022-2025.")
    print("  That gap needs an explicit, documented decision at Stage E, not silent interpolation.")


if __name__ == "__main__":
    main()
