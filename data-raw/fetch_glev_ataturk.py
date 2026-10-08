"""
data-raw/fetch_glev_ataturk.py

Stage C.5 (Ataturk_Master_Checklist_v2.md): independent benchmark from the
Global Lake Evaporation Volume (GLEV) dataset, Zenodo
10.5281/zenodo.4646621 (Zhao et al. 2022, Nature Communications).

What was verified before writing this (7 Oct 2026, from the Zenodo record
and the paper, not from memory):
  - Three CSVs, ~4 GB each, CC BY 4.0: 0_evaporation_rate.csv (mm/day),
    1_openwater_area.csv (m2, lake ice already removed),
    2_evaporation_volume.csv (thousand m3 per month).
  - Each has 1,427,687 rows x 409 columns: an identifier, then 408 monthly
    values, January 1985 to December 2018. Rows are keyed by the
    HydroLAKES v1.0 `Hylak_id`; the files carry NO coordinates or names,
    so Ataturk has to be found in HydroLAKES first.
  - Method (paper): Penman combination equation with lake heat storage,
    forcing from TerraClimate/ERA5/GLDAS, monthly Landsat-based area
    (GSW) inside HydroLAKES polygons; volume = rate x area x (1 - ice).
    Total uncertainty stated as 9.93%.
  - HydroLAKES attributes (names confirmed): Hylak_id, Lake_name, Country,
    Lake_type (2 = reservoir), Grand_id, Lake_area (km2), Vol_total (mcm),
    Elevation (m), Pour_long, Pour_lat. Pour-point layer is a separate
    79 MB shapefile zip at data.hydrosheds.org.

NOT verified, so the script checks it at run time instead of assuming:
  - whether the CSVs have a header row (detected from the first line),
  - that the 408 value columns are in chronological order Jan 1985 ->
    Dec 2018 (taken from the record's description; the printed column
    names are shown so you can eyeball them),
  - that the three files agree with each other: volume must equal
    rate x days-in-month x area for every month (this also confirms the
    units read from the record). The script stops if they do not.
  - which HydroLAKES row is Ataturk: found two independent ways (name
    match, and a pour point inside the JRC polygon's bounding box with a
    large area). If they disagree it stops and prints the candidates.

USE IN THE PAPER: GLEV is an independent benchmark, like NASA POWER is an
independent forcing member. It is not a calibration target.

PREREQUISITES: pip install pyshp shapely   (shapely is already needed by
fetch_modis_lst.py). Everything else is the standard library.
DISK/NETWORK: up to ~12 GB streamed from Zenodo and ~80 MB for
HydroLAKES. The big CSVs are streamed and not stored; the script stops
reading each file as soon as it finds the row. If you already downloaded
them, put them in data-raw/glev/ and they are read from there instead.
"""

import calendar
import csv
import urllib.error
import json
import os
import sys
import time
import urllib.request
import zipfile

import shapefile  # pyshp
from shapely.geometry import shape

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(SCRIPT_DIR, "..", "data")
HL_DIR = os.path.join(SCRIPT_DIR, "hydrolakes")   # git-ignored (data-raw/*)
GLEV_DIR = os.path.join(SCRIPT_DIR, "glev")       # git-ignored
os.makedirs(DATA_DIR, exist_ok=True)
os.makedirs(HL_DIR, exist_ok=True)
os.makedirs(GLEV_DIR, exist_ok=True)

POLYGON_GEOJSON_PATH = os.path.join(DATA_DIR, "ataturk_jrc_max_extent_polygon.geojson")
OUT_CSV = os.path.join(DATA_DIR, "glev_ataturk_1985-2018.csv")

HL_URL = "https://data.hydrosheds.org/file/hydrolakes/HydroLAKES_points_v10_shp.zip"
ZENODO = "https://zenodo.org/records/4646621/files/{name}?download=1"
FILES = {
    "rate":   "0_evaporation_rate.csv",
    "area":   "1_openwater_area.csv",
    "volume": "2_evaporation_volume.csv",
}
N_MONTHS = 408  # Jan 1985 .. Dec 2018

# data.hydrosheds.org returned HTTP 403 to Python's default urllib
# User-Agent on the first real run (8 Oct 2026). Sending a browser-style
# header is the usual fix; if it still 403s, download the zip in a
# browser and drop it in data-raw/hydrolakes/ - the script skips the
# download when the file is already there.
UA = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 "
                    "(KHTML, like Gecko) Chrome/124.0 Safari/537.36"}


def http_open(url):
    return urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=120)
BBOX_MARGIN_DEG = 0.15
MIN_AREA_KM2 = 300.0   # Ataturk is the only lake this large near the dam


# ---------------------------------------------------------------
# 1. Find Ataturk's Hylak_id in HydroLAKES
# ---------------------------------------------------------------
def download(url, dest):
    print(f"  downloading {url}")
    with http_open(url) as r, open(dest, "wb") as f:
        total = 0
        while True:
            chunk = r.read(1 << 20)
            if not chunk:
                break
            f.write(chunk)
            total += len(chunk)
    print(f"  saved {total / 1e6:.1f} MB -> {dest}")


def find_hylak_id():
    zpath = os.path.join(HL_DIR, "HydroLAKES_points_v10_shp.zip")
    if not os.path.exists(zpath):
        try:
            download(HL_URL, zpath)
        except urllib.error.HTTPError as e:
            if os.path.exists(zpath):
                os.remove(zpath)
            sys.exit(f"HydroLAKES download refused (HTTP {e.code}). Download "
                     f"HydroLAKES_points_v10_shp.zip from "
                     f"https://www.hydrosheds.org/hydrolakes in a browser, save it as\n  {zpath}\n"
                     f"and run this script again.")
    with zipfile.ZipFile(zpath) as z:
        dbfs = [n for n in z.namelist() if n.lower().endswith(".dbf")]
        if not dbfs:
            sys.exit("No .dbf found inside the HydroLAKES zip.")
        z.extractall(HL_DIR)
    dbf_path = os.path.join(HL_DIR, dbfs[0])

    with open(dbf_path, "rb") as f:
        sf = shapefile.Reader(dbf=f, encoding="utf-8", encodingErrors="replace")
        names = [fld[0] for fld in sf.fields[1:]]
        print(f"  HydroLAKES fields: {names}")
        needed = {"Hylak_id", "Lake_name", "Lake_area", "Pour_long", "Pour_lat"}
        if not needed <= set(names):
            sys.exit(f"Expected fields missing: {sorted(needed - set(names))}")
        recs = [r.as_dict() for r in sf.records()]
    print(f"  {len(recs):,} HydroLAKES records")

    with open(POLYGON_GEOJSON_PATH, encoding="utf-8") as f:
        gj = json.load(f)
    geom = gj["geometry"] if gj.get("type") == "Feature" else gj
    minx, miny, maxx, maxy = shape(geom).bounds
    minx -= BBOX_MARGIN_DEG; miny -= BBOX_MARGIN_DEG
    maxx += BBOX_MARGIN_DEG; maxy += BBOX_MARGIN_DEG

    def fnum(x):
        try:
            return float(x)
        except (TypeError, ValueError):
            return float("nan")

    by_name = [r for r in recs if "atat" in str(r["Lake_name"]).lower()]
    by_place = [r for r in recs
                if minx <= fnum(r["Pour_long"]) <= maxx
                and miny <= fnum(r["Pour_lat"]) <= maxy
                and fnum(r["Lake_area"]) >= MIN_AREA_KM2]

    def show(title, rows):
        print(f"  {title}: {len(rows)} row(s)")
        for r in rows:
            print("    Hylak_id=%s name=%r country=%r area=%.1f km2 vol=%s mcm "
                  "elev=%s m pour=(%.4f, %.4f) Grand_id=%s Lake_type=%s" % (
                      int(fnum(r["Hylak_id"])), r.get("Lake_name"), r.get("Country"),
                      fnum(r["Lake_area"]), r.get("Vol_total"), r.get("Elevation"),
                      fnum(r["Pour_long"]), fnum(r["Pour_lat"]),
                      r.get("Grand_id"), r.get("Lake_type")))

    show("name contains 'Atat'", by_name)
    show(f"pour point inside JRC bbox (+{BBOX_MARGIN_DEG} deg), area >= {MIN_AREA_KM2:.0f} km2", by_place)

    ids_name = {int(fnum(r["Hylak_id"])) for r in by_name}
    ids_place = {int(fnum(r["Hylak_id"])) for r in by_place}
    common = ids_name & ids_place
    if len(common) == 1:
        hid = common.pop()
    elif len(ids_place) == 1 and not ids_name:
        hid = next(iter(ids_place))
        print("  NOTE: no row has 'Atat' in its name; using the single large lake "
              "whose pour point is at the dam. Check the attributes above.")
    else:
        sys.exit("Could not identify Ataturk uniquely - see the candidates above. "
                 "Set HYLAK_ID_OVERRIDE below once you have chosen one.")
    print(f"  -> using Hylak_id {hid}")
    return hid


HYLAK_ID_OVERRIDE = None  # set to an int to skip the search


# ---------------------------------------------------------------
# 2. Pull that row out of each GLEV CSV (streamed, early exit)
# ---------------------------------------------------------------
def open_source(fname):
    local = os.path.join(GLEV_DIR, fname)
    if os.path.exists(local):
        print(f"  reading local copy {local}")
        return open(local, "rb")
    url = ZENODO.format(name=fname)
    print(f"  streaming {url}")
    return http_open(url)


def first_field_int(raw):
    head = raw.split(b",", 1)[0].strip().strip(b'"')
    try:
        return int(float(head))
    except ValueError:
        return None


def extract_row(fname, hid, attempts=3):
    for attempt in range(1, attempts + 1):
        try:
            src = open_source(fname)
            header, found, nread, t0 = None, None, 0, time.time()
            with src:
                for i, raw in enumerate(src):
                    nread += len(raw)
                    if i == 0:
                        if first_field_int(raw) is None:
                            header = raw.decode("utf-8", "replace").rstrip("\r\n").split(",")
                            continue
                    if first_field_int(raw) == hid:
                        found = raw.decode("utf-8", "replace").rstrip("\r\n").split(",")
                        break
                    if i % 100000 == 0 and i:
                        print(f"    ...{i:,} rows, {nread / 1e9:.2f} GB, "
                              f"{(time.time() - t0) / 60:.1f} min")
            if found is None:
                sys.exit(f"Hylak_id {hid} not found in {fname}.")
            return header, found
        except (OSError, urllib.error.URLError) as e:
            print(f"    attempt {attempt}/{attempts} failed: {e}")
            if attempt == attempts:
                raise
            time.sleep(10)


def to_floats(fields):
    out = []
    for x in fields:
        x = x.strip().strip('"')
        try:
            out.append(float(x))
        except ValueError:
            out.append(float("nan"))
    return out


def main():
    t_start = time.time()
    print("Step 1: locate Ataturk in HydroLAKES")
    hid = HYLAK_ID_OVERRIDE or find_hylak_id()

    series = {}
    for key, fname in FILES.items():
        print(f"\nStep 2: {fname}")
        header, row = extract_row(fname, hid)
        if len(row) != N_MONTHS + 1:
            sys.exit(f"{fname}: expected {N_MONTHS + 1} columns, got {len(row)}.")
        if header:
            print(f"    header present; first columns: {header[:3]} ... last: {header[-2:]}")
        else:
            print("    no header row detected")
        series[key] = to_floats(row[1:])
        print(f"    row found: {len(series[key])} monthly values")

    # --- build the monthly table and cross-check the three files ----
    rows, worst = [], 0.0
    for i in range(N_MONTHS):
        year, month = 1985 + i // 12, i % 12 + 1
        days = calendar.monthrange(year, month)[1]
        rate, area_m2, vol_k = series["rate"][i], series["area"][i], series["volume"][i]
        if all(v == v for v in (rate, area_m2, vol_k)) and vol_k > 0:
            expect = rate * days * area_m2 / 1e6   # thousand m3/month
            worst = max(worst, abs(expect - vol_k) / vol_k)
        rows.append({
            "date": f"{year}-{month:02d}-01",
            "hylak_id": hid,
            "evap_rate_mm_day": rate,
            "openwater_area_km2": area_m2 / 1e6,
            "evap_volume_mcm_month": vol_k / 1e3,
        })
    print(f"\nCross-check: volume vs rate x days x area, worst relative difference "
          f"= {100 * worst:.3f}%")
    if worst > 0.02:
        sys.exit("The three files do not agree to within 2% - the column order or "
                 "units assumed above are wrong. Not writing output.")

    with open(OUT_CSV, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"Wrote {len(rows)} rows to {OUT_CSV}")

    # --- summary to eyeball ---------------------------------------
    def mean(xs):
        xs = [x for x in xs if x == x]
        return sum(xs) / len(xs) if xs else float("nan")

    yearly_depth, yearly_vol = {}, {}
    for r in rows:
        y = int(r["date"][:4]); days = calendar.monthrange(y, int(r["date"][5:7]))[1]
        yearly_depth.setdefault(y, 0.0)
        yearly_depth[y] += r["evap_rate_mm_day"] * days
        yearly_vol.setdefault(y, 0.0)
        yearly_vol[y] += r["evap_volume_mcm_month"]
    print("\n=== GLEV summary for this reservoir (1985-2018) ===")
    print(f"  mean open-water area : {mean([r['openwater_area_km2'] for r in rows]):8.1f} km2")
    print(f"  mean evaporation     : {mean(list(yearly_depth.values())):8.0f} mm/yr")
    print(f"  mean annual volume   : {mean(list(yearly_vol.values())):8.0f} million m3/yr "
          f"({mean(list(yearly_vol.values())) / 1000:.2f} km3/yr)")
    for y in (1985, 1995, 2005, 2015, 2018):
        print(f"    {y}: {yearly_depth[y]:6.0f} mm, {yearly_vol[y]:7.0f} million m3")
    print(f"\nTotal run time: {(time.time() - t_start) / 60:.1f} minutes")


if __name__ == "__main__":
    main()
