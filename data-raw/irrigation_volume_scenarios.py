"""Annual irrigation-withdrawal SCENARIOS from satellite summer-green area x duty.

Reads  data/irrigated_area_summer_1993-2025.csv   (fetch_irrigated_area.py)
Writes data/irrigation_withdrawal_scenarios_1993-2025.csv

This is a scale estimate of WITHDRAWAL, not net consumption, and the area is
'July-August green (NDVI >= 0.40) area inside three rectangles', NOT the
DSI command-area. See Ataturk_Master_Checklist_v2.md, C.6, for the caveats.

  low  = Harran box area, capped at the district size implied by Ozdogan et al.
         (2006): 101,524 ha = ~76% of district land -> ~133,584 ha;
         x DUTY_LOW.
  high = Harran + Suruc + Ceylanpinar-Mardin boxes, uncapped; x DUTY_HIGH.
Duty 10,000-12,000 m3/ha/yr: Goor, Alia, van der Zaag & Tilmant (2007), IAHS 315.
"""
import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "data" / "irrigated_area_summer_1993-2025.csv"
OUT = ROOT / "data" / "irrigation_withdrawal_scenarios_1993-2025.csv"

DUTY_LOW, DUTY_HIGH = 10_000.0, 12_000.0           # m3/ha/yr
HARRAN_DISTRICT_HA = 101_524 / 0.76                # derived, ~133,584 ha
BOXES = ["HARRAN", "SURUC", "CEYLANPINAR_MARDIN"]
M3_PER_KM3 = 1e9

rows = {}
with open(SRC, newline="") as f:
    for r in csv.DictReader(f):
        if r["flag"] not in ("", "NA", "None"):
            raise SystemExit(f"flagged row {r['year']} {r['box']}: {r['flag']}")
        rows[(r["box"], int(r["year"]))] = float(r["area_est_ha_T40"])

years = sorted({y for (_, y) in rows})
assert all((b, y) in rows for b in BOXES for y in years), "missing box-year"

with open(OUT, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["year", "harran_ha", "suruc_ha", "ceylanpinar_mardin_ha", "three_box_ha",
                "harran_capped_ha", "low_km3", "high_km3"])
    for y in years:
        h, s, c = (rows[(b, y)] for b in BOXES)
        hc = min(h, HARRAN_DISTRICT_HA)
        w.writerow([y, round(h), round(s), round(c), round(h + s + c), round(hc),
                    round(hc * DUTY_LOW / M3_PER_KM3, 3),
                    round((h + s + c) * DUTY_HIGH / M3_PER_KM3, 3)])
print(f"wrote {OUT} ({len(years)} years); Harran cap = {HARRAN_DISTRICT_HA:,.0f} ha")
