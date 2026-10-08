"""
data-raw/fetch_era5land_wind_at_stations.py   (8 Oct 2026)

Wind test, step 2 (Ataturk_Master_Checklist_v2.md, C.4 wind block). The
station check showed observed regional wind (~2.6 m/s at the three
homogeneous stations) lies between ERA5-Land (reservoir mean 1.99) and
POWER (3.53). That compared a RESERVOIR-mean ERA5-Land value with
STATION observations, so lake and land were mixed. This script extracts
ERA5-Land 10 m wind AT THE STATION COORDINATES so the comparison is
like-for-like.

Definition matches the project's own wspd10: sqrt(u10^2 + v10^2) at every
HOUR first, then the mean over the UTC day (NOT the magnitude of the
daily-mean components, which was the 29 Sep bug).

Verified (Earth Engine catalog page, 8 Oct 2026): collection
ECMWF/ERA5_LAND/HOURLY; bands u_component_of_wind_10m and
v_component_of_wind_10m, m/s; pixel size 11132 m; coverage into Oct 2026.
The script re-checks band names at run time and stops if they differ.

Points: the seven ISD stations (coordinates as printed by
wind_station_check.R from isd-history.csv) plus the reservoir centroid
(sanity: its long-run mean should be in the neighbourhood of the project's
reservoir-mean 1.99 m/s; it is one cell, not the lake-weighted mean).

Output: data/era5land_wind_at_stations_1993-2025.csv
        columns: date, point, name, lon, lat, wspd10_era5l, n_hours
Run:    python data-raw/fetch_era5land_wind_at_stations.py
Needs:  earthengine-api (already installed for the MODIS/JRC fetches).
Takes:  roughly 10-30 minutes (33 yearly requests, retried on timeout).
"""
import calendar
import csv
import os
import sys
import time

import ee

EE_PROJECT = "ataturk-evaporation"
COLLECTION = "ECMWF/ERA5_LAND/HOURLY"
U_BAND = "u_component_of_wind_10m"
V_BAND = "v_component_of_wind_10m"
YEAR_START, YEAR_END = 1993, 2025          # inclusive
GETINFO_MAX_ATTEMPTS = 3
GETINFO_RETRY_WAIT_S = 15

POINTS = [  # id, name, lon, lat
    ("172700-99999", "SANLIURFA",     38.767, 37.133),
    ("171970-99999", "TULGA",         38.254, 38.354),
    ("172000-99999", "ERHAC/MALATYA", 38.091, 38.435),
    ("172600-99999", "OGUZELI/GAZIANTEP", 37.479, 36.947),
    ("172020-99999", "ELAZIG",        39.291, 38.607),
    ("172800-99999", "DIYARBAKIR",    40.201, 37.894),
    ("172550-99999", "KAHRAMANMARAS", 36.933, 37.600),
    ("RESERVOIR_CENTROID", "reservoir centroid", 38.5878, 37.6032),
]

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "data", "era5land_wind_at_stations_1993-2025.csv")


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


def main():
    t0 = time.time()
    ee.Authenticate()
    ee.Initialize(project=EE_PROJECT)

    era = ee.ImageCollection(COLLECTION)
    bands = ee.Image(era.first()).bandNames().getInfo()
    for b in (U_BAND, V_BAND):
        if b not in bands:
            sys.exit(f"Band {b!r} not found in {COLLECTION}. Bands: {bands[:12]}...")
    era = era.select([U_BAND, V_BAND])
    proj = ee.Image(era.first()).select(U_BAND).projection()
    print("ERA5-Land native pixel (m):", proj.nominalScale().getInfo())

    pts = ee.FeatureCollection([
        ee.Feature(ee.Geometry.Point([lon, lat]), {"point": pid, "name": name, "lon": lon, "lat": lat})
        for pid, name, lon, lat in POINTS])

    rows = []
    for year in range(YEAR_START, YEAR_END + 1):
        n_days = 366 if calendar.isleap(year) else 365
        start = ee.Date(f"{year}-01-01")

        def per_day(i):
            d0 = start.advance(ee.Number(i), "day")
            hourly = era.filterDate(d0, d0.advance(1, "day"))
            speed = hourly.map(lambda im: im.select(U_BAND).hypot(im.select(V_BAND)).rename("ws"))
            return (speed.mean().set("date", d0.format("YYYY-MM-dd"))
                    .set("n_hours", hourly.size()))

        daily = ee.ImageCollection.fromImages(ee.List.sequence(0, n_days - 1).map(per_day))

        def sample(img):
            fc = img.reduceRegions(collection=pts, reducer=ee.Reducer.mean(),
                                   crs=proj, scale=11132)
            return fc.map(lambda f: f.set("date", img.get("date")).set("n_hours", img.get("n_hours")))

        info = getinfo_with_retry(daily.map(sample).flatten(), f"{year}")
        feats = info["features"]
        if year == YEAR_START:
            print("  property check (first feature):", sorted(feats[0]["properties"].keys()))
        for f in feats:
            pr = f["properties"]
            if "mean" not in pr or pr["mean"] is None:
                sys.exit(f"{year}: no 'mean' for {pr.get('point')} on {pr.get('date')}; properties {pr}")
            rows.append((pr["date"], pr["point"], pr["name"], pr["lon"], pr["lat"], pr["mean"], int(pr["n_hours"])))
        short = sum(1 for f in feats if int(f["properties"]["n_hours"]) != 24)
        print(f"  {year}: {len(feats)} point-days"
              + (f"  ({short} with n_hours != 24)" if short else ""))

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["date", "point", "name", "lon", "lat", "wspd10_era5l", "n_hours"])
        w.writerows(rows)
    print(f"\nWrote {len(rows)} rows to {os.path.abspath(OUT)}")

    print("\nLong-run mean of daily-mean hourly speed (m/s):")
    by = {}
    for r in rows:
        by.setdefault(r[1], []).append(r[5])
    for pid, name, _, _ in POINTS:
        v = by[pid]
        print(f"  {pid:22s} {name:20s} n={len(v)}  mean={sum(v) / len(v):.2f}  max={max(v):.1f}")
    print("  (project reservoir-mean wspd10, lake-weighted: 1.99)")
    print(f"\nTotal run time: {(time.time() - t0) / 60:.1f} minutes")


if __name__ == "__main__":
    main()
