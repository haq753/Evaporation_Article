"""
data-raw/fetch_irrigated_area.py   (9 Oct 2026)

Stage C.6 follow-up (Ataturk_Master_Checklist_v2.md): an annual irrigated-area
series for the three schemes fed by the Sanliurfa tunnels, 1993-2025, to extend
Ozdogan et al. (2006, Water Resources Management 20:467-488), whose Harran
Plain record ends in 2002. Used to turn the checklist's dated irrigation range
(duty x area) into an annual line.

METHOD (deliberately simple and checkable; same logic as Ozdogan's "summer
irrigated crops"): summer in this region is essentially rain-free, so land that
is GREEN in July-August is irrigated (or riparian/wetland/irrigated trees).
For each year: take every clear Landsat surface-reflectance observation from
1 Jul to 31 Aug, compute NDVI, take the per-pixel MAXIMUM, and count the area
with max NDVI >= T (primary T = 0.40; 0.35/0.45/0.50 reported as a sensitivity
band; T was fixed a priori, not tuned).

VERIFIED against the Earth Engine catalog pages on 9 Oct 2026 (not recalled):
  LANDSAT/LT05/C02/T1_L2  1984-03-16 to 2012-05-05  red SR_B3  NIR SR_B4
  LANDSAT/LE07/C02/T1_L2  1999-05-28 to 2024-01-19  red SR_B3  NIR SR_B4
  LANDSAT/LC08/C02/T1_L2  2013-03-18 onwards         red SR_B4  NIR SR_B5
  LANDSAT/LC09/C02/T1_L2  2021-10-31 onwards         red SR_B4  NIR SR_B5
  SR scale 2.75e-05, offset -0.2 (all four).
  QA_PIXEL bits: 0 fill, 1 dilated cloud, 2 cirrus (LC08/LC09 only), 3 cloud,
  4 cloud shadow, 5 snow, 6 clear (inverted sense), 7 water.
  Consequence: 1 Jul-31 Aug 2012 has ONLY Landsat 7 (SLC-off striping), and 2013
  onward has two or more sensors.

WHAT THIS DOES NOT DO (state in Methods):
  * It measures GREEN summer area, not water source. Groundwater-irrigated land
    is counted too, and in the Harran Plain wells are used alongside the canals.
  * Orchards, riparian strips and gardens that are green in summer are counted.
  * The three boxes below are PROVISIONAL rectangles (I do not have the DSI
    command-area polygons). Check them on the preview PNGs before using the
    totals. Only HARRAN can be checked against Ozdogan (2002: 101,524 ha).
  * Landsat 5 coverage over Turkey in the 1990s is patchy; years with little
    clear data are flagged, not filled. area_est_ha_T40 = green area on valid
    pixels / valid fraction (assumes the valid part is representative).
  * Landsat 5/7 and 8/9 NDVI differ slightly; per-sensor areas are written to a
    second file so the step at the 2013 sensor change can be inspected.

OUTPUTS
  data/irrigated_area_summer_1993-2025.csv          year x box
  data/irrigated_area_by_sensor_1993-2025.csv       year x box x sensor (T = 0.40)
  data/irrigation_preview_<year>.png                max-NDVI map with the boxes
The script RESUMES: rows already in the first file are skipped on re-run.
Run:  python data-raw/fetch_irrigated_area.py       (Earth Engine, ~30-90 min)
"""
import csv
import os
import sys
import time
import urllib.request
from datetime import datetime, timezone

import ee

EE_PROJECT = "ataturk-evaporation"
YEAR_START, YEAR_END = 1993, 2025
WINDOW_START_MD, WINDOW_END_MD = "07-01", "09-01"      # end date is exclusive
THRESHOLDS = [0.35, 0.40, 0.45, 0.50]
MAIN_T = 0.40
SCALE_M = 30
CRS = "EPSG:32637"          # UTM 37N covers 36-42 E, i.e. all three boxes
LOW_COVERAGE_FRAC = 0.80
PREVIEW_YEARS = [2002, 2010, 2021]
GETINFO_MAX_ATTEMPTS = 3
GETINFO_RETRY_WAIT_S = 20

BOXES = {   # name: (west, south, east, north) in degrees; PROVISIONAL
    "HARRAN":             (38.75, 36.70, 39.60, 37.20),
    "SURUC":              (38.25, 36.80, 38.70, 37.20),
    "CEYLANPINAR_MARDIN": (39.85, 36.65, 40.90, 37.35),
}

SENSORS = {
    "LT05": dict(coll="LANDSAT/LT05/C02/T1_L2", red="SR_B3", nir="SR_B4", cirrus=False),
    "LE07": dict(coll="LANDSAT/LE07/C02/T1_L2", red="SR_B3", nir="SR_B4", cirrus=False),
    "LC08": dict(coll="LANDSAT/LC08/C02/T1_L2", red="SR_B4", nir="SR_B5", cirrus=True),
    "LC09": dict(coll="LANDSAT/LC09/C02/T1_L2", red="SR_B4", nir="SR_B5", cirrus=True),
}

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "..", "data")
OUT_MAIN = os.path.join(DATA, "irrigated_area_summer_1993-2025.csv")
OUT_SENS = os.path.join(DATA, "irrigated_area_by_sensor_1993-2025.csv")

MAIN_FIELDS = (["year", "box", "total_ha", "valid_ha", "valid_frac", "n_dates", "sensors"]
               + [f"green_ha_T{int(round(t * 100))}" for t in THRESHOLDS]
               + ["area_est_ha_T40", "flag"])
SENS_FIELDS = ["year", "box", "sensor", "n_images", "valid_frac", "green_ha_T40"]


# ---------------------------------------------------------------- pure helpers
def tkey(t):
    return f"T{int(round(t * 100))}"


def postprocess(year, box, props, sensors_present, n_dates):
    """Turn the summed band properties of one (year, box) into output rows.
    Pure function: tested without Earth Engine."""
    total = float(props["total_ha"])
    valid = float(props.get("all_valid") or 0.0)
    vfrac = valid / total if total > 0 else 0.0
    greens = {tkey(t): float(props.get(f"all_g{tkey(t)}") or 0.0) for t in THRESHOLDS}
    if vfrac <= 0:
        est, flag = "", "NO_DATA"
    else:
        est = greens["T40"] / vfrac
        flag = "LOW_COVERAGE" if vfrac < LOW_COVERAGE_FRAC else ""
    row = {"year": year, "box": box, "total_ha": round(total, 1), "valid_ha": round(valid, 1),
           "valid_frac": round(vfrac, 4), "n_dates": n_dates, "sensors": "+".join(sensors_present),
           "area_est_ha_T40": "" if est == "" else round(est, 1), "flag": flag}
    for t in THRESHOLDS:
        row[f"green_ha_{tkey(t)}"] = round(greens[tkey(t)], 1)
    srows = []
    for s in sensors_present:
        sv = float(props.get(f"{s}_valid") or 0.0)
        srows.append({"year": year, "box": box, "sensor": s,
                      "n_images": "", "valid_frac": round(sv / total, 4) if total > 0 else 0.0,
                      "green_ha_T40": round(float(props.get(f"{s}_gT40") or 0.0), 1)})
    return row, srows


def getinfo_with_retry(obj, label):
    last = None
    for attempt in range(1, GETINFO_MAX_ATTEMPTS + 1):
        try:
            return obj.getInfo()
        except ee.ee_exception.EEException as e:
            last = e
            if "timed out" in str(e).lower() and attempt < GETINFO_MAX_ATTEMPTS:
                print(f"    {label}: timed out (attempt {attempt}/{GETINFO_MAX_ATTEMPTS}); "
                      f"waiting {GETINFO_RETRY_WAIT_S}s and retrying.")
                time.sleep(GETINFO_RETRY_WAIT_S)
                continue
            raise
    raise last


# ---------------------------------------------------------------- EE building blocks
def ndvi_collection(sensor, year, geom):
    cfg = SENSORS[sensor]
    col = (ee.ImageCollection(cfg["coll"]).filterBounds(geom)
           .filterDate(f"{year}-{WINDOW_START_MD}", f"{year}-{WINDOW_END_MD}"))

    def to_ndvi(img):
        qa = img.select("QA_PIXEL")
        bad = (qa.bitwiseAnd(1 << 0).neq(0)          # fill
               .Or(qa.bitwiseAnd(1 << 1).neq(0))     # dilated cloud
               .Or(qa.bitwiseAnd(1 << 3).neq(0))     # cloud
               .Or(qa.bitwiseAnd(1 << 4).neq(0)))    # cloud shadow
        if cfg["cirrus"]:
            bad = bad.Or(qa.bitwiseAnd(1 << 2).neq(0))
        red = img.select(cfg["red"]).multiply(2.75e-05).add(-0.2)
        nir = img.select(cfg["nir"]).multiply(2.75e-05).add(-0.2)
        ndvi = nir.subtract(red).divide(nir.add(red)).rename("ndvi")
        ok = red.gt(0).And(nir.gt(0))                # physically sensible reflectance
        return ndvi.updateMask(bad.Not()).updateMask(ok)
    return col, col.map(to_ndvi)


def area_bands(ndvi, prefix, thresholds):
    px = ee.Image.pixelArea().divide(10000.0)        # hectares
    bands = [px.updateMask(ndvi.mask()).rename(f"{prefix}_valid")]
    for t in thresholds:
        bands.append(px.updateMask(ndvi.gte(t)).rename(f"{prefix}_g{tkey(t)}"))
    return bands


def box_geom(b):
    return ee.Geometry.Rectangle(list(b), proj="EPSG:4326", geodesic=False)


def main():
    t0 = time.time()
    ee.Authenticate()
    ee.Initialize(project=EE_PROJECT)
    os.makedirs(DATA, exist_ok=True)

    union = ee.Geometry.Rectangle([min(b[0] for b in BOXES.values()), min(b[1] for b in BOXES.values()),
                                   max(b[2] for b in BOXES.values()), max(b[3] for b in BOXES.values())],
                                  proj="EPSG:4326", geodesic=False)
    box_geoms = {k: box_geom(v) for k, v in BOXES.items()}

    done = set()
    if os.path.exists(OUT_MAIN):
        with open(OUT_MAIN, newline="") as fh:
            for r in csv.DictReader(fh):
                done.add((int(r["year"]), r["box"]))
    new_main = not os.path.exists(OUT_MAIN)
    new_sens = not os.path.exists(OUT_SENS)
    fm = open(OUT_MAIN, "a", newline="")
    fs = open(OUT_SENS, "a", newline="")
    wm = csv.DictWriter(fm, fieldnames=MAIN_FIELDS)
    ws = csv.DictWriter(fs, fieldnames=SENS_FIELDS)
    if new_main:
        wm.writeheader()
    if new_sens:
        ws.writeheader()

    for year in range(YEAR_START, YEAR_END + 1):
        todo = [b for b in BOXES if (year, b) not in done]
        # --- which sensors / dates exist this summer (metadata only, cheap)
        cols, mapped, times = {}, {}, {}
        for s in SENSORS:
            cols[s], mapped[s] = ndvi_collection(s, year, union)
            times[s] = cols[s].aggregate_array("system:time_start")
        meta = getinfo_with_retry(ee.Dictionary({s: times[s] for s in SENSORS}), f"{year} metadata")
        present = [s for s in SENSORS if meta[s]]
        dates = sorted({datetime.fromtimestamp(ms / 1000, tz=timezone.utc).strftime("%Y-%m-%d")
                        for s in present for ms in meta[s]})
        if not present:
            print(f"{year}: no Landsat scenes in the window")
        if not todo and year not in PREVIEW_YEARS:
            print(f"{year}: already done")
            continue

        if present:
            merged = None
            for s in present:
                merged = mapped[s] if merged is None else merged.merge(mapped[s])
            all_ndvi = merged.max().rename("ndvi")
            px = ee.Image.pixelArea().divide(10000.0).rename("total_ha")
            bands = [px] + area_bands(all_ndvi, "all", THRESHOLDS)
            for s in present:
                bands += area_bands(mapped[s].max().rename("ndvi"), s, [MAIN_T])
            img = ee.Image.cat(bands)

            for b in todo:
                fc = img.reduceRegions(collection=ee.FeatureCollection([ee.Feature(box_geoms[b], {"box": b})]),
                                       reducer=ee.Reducer.sum(), scale=SCALE_M, crs=CRS, tileScale=4)
                feats = getinfo_with_retry(fc, f"{year} {b}")["features"]
                props = feats[0]["properties"]
                row, srows = postprocess(year, b, props, present, len(dates))
                for sr in srows:
                    sr["n_images"] = len(meta[sr["sensor"]])
                wm.writerow(row)
                for sr in srows:
                    ws.writerow(sr)
                fm.flush(); fs.flush()
                print(f"{year} {b:19s} sensors {row['sensors']:<18s} dates {row['n_dates']:>2d}  "
                      f"valid {row['valid_frac']:.2f}  green(T40) {row['green_ha_T40']:>10,.0f} ha  "
                      f"est {row['area_est_ha_T40'] if row['area_est_ha_T40'] == '' else format(row['area_est_ha_T40'], ',.0f')}"
                      f"  {row['flag']}")
        else:
            for b in todo:
                total = float(getinfo_with_retry(
                    ee.Image.pixelArea().divide(10000.0).reduceRegion(
                        ee.Reducer.sum(), box_geoms[b], SCALE_M, crs=CRS, maxPixels=int(1e10)), f"{year} {b} area")["area"])
                row, _ = postprocess(year, b, {"total_ha": total}, [], 0)
                wm.writerow(row); fm.flush()
                print(f"{year} {b:19s} NO_DATA")

        if year in PREVIEW_YEARS and present:
            try:
                base = ee.Image.constant(0).visualize(min=0, max=1, palette=["808080"])
                veg = all_ndvi.visualize(min=0.0, max=0.8, palette=["ffffcc", "c2e699", "78c679", "31a354", "006837"])
                outlines = ee.FeatureCollection([ee.Feature(g) for g in box_geoms.values()]).style(
                    color="ff0000", fillColor="00000000", width=2)
                url = base.blend(veg).blend(outlines).getThumbURL(
                    {"region": union, "dimensions": 1400, "format": "png", "crs": "EPSG:4326"})
                path = os.path.join(DATA, f"irrigation_preview_{year}.png")
                urllib.request.urlretrieve(url, path)
                print(f"    preview saved: {os.path.abspath(path)}")
            except Exception as e:  # preview is a convenience, never fatal
                print(f"    preview for {year} failed: {e}")

    fm.close(); fs.close()
    print(f"\nWrote {os.path.abspath(OUT_MAIN)} and {os.path.abspath(OUT_SENS)}")
    print(f"Total run time: {(time.time() - t0) / 60:.1f} minutes")


if __name__ == "__main__":
    main()
