"""
Per-metric percentile grid over the PROTAC-DB reference set.

Powers UI context like "MW is at the 82nd percentile of known PROTACs".

Kept DELIBERATELY SEPARATE from calibration.json (which drives scoring): building
this never perturbs the scoring thresholds, so webapp is unaffected. Additive file,
additive module — nothing in the scoring path imports it.

Build:
  <venv>/python -m engine.distribution data/protac_smiles.smi --source "PROTAC-DB"
"""
import argparse
import json
import os

import numpy as np
from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")
from .metrics import METRICS

_PATH = os.path.join(os.path.dirname(__file__), "distribution.json")
_GRID = None
_META = None


def build(smiles_file, source, out=_PATH):
    from .calibrate import read_smiles
    smis = read_smiles(smiles_file)
    cols = {k: [] for k in METRICS}
    n_ok = 0
    for smi in smis:
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        n_ok += 1
        for k, (fn, _b, _w, _h) in METRICS.items():
            try:
                v = fn(m)
            except Exception:
                v = None
            if v is not None:
                cols[k].append(float(v))

    grid = {}
    for k, vals in cols.items():
        if len(vals) < 50:
            continue
        a = np.array(vals, dtype=float)
        grid[k] = [round(float(np.percentile(a, p)), 3) for p in range(101)]  # anchors at p=0..100

    obj = {"meta": {"source": source, "n_valid": n_ok}, "grid": grid}
    with open(out, "w") as fh:
        json.dump(obj, fh)
    return obj


def _load():
    global _GRID, _META
    if _GRID is not None:
        return
    if not os.path.exists(_PATH):
        _GRID = {}
        return
    with open(_PATH) as fh:
        d = json.load(fh)
    _GRID = d.get("grid", {})
    _META = d.get("meta")


def percentile_of(metric, value):
    """Approx percentile (0-100) of `value` within the reference set, or None."""
    _load()
    if value is None or metric not in _GRID:
        return None
    anchors = _GRID[metric]                 # values at percentiles 0..100 (non-decreasing)
    if value <= anchors[0]:
        return 0
    if value >= anchors[100]:
        return 100
    for p in range(100):
        if anchors[p] <= value <= anchors[p + 1]:
            span = anchors[p + 1] - anchors[p]
            frac = 0.0 if span == 0 else (value - anchors[p]) / span
            return int(round(p + frac))
    return 100


def meta():
    _load()
    return _META


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("smiles_file")
    ap.add_argument("--source", default="PROTAC-DB")
    a = ap.parse_args()
    obj = build(a.smiles_file, a.source)
    print(f"built distribution grid: n_valid={obj['meta']['n_valid']}, metrics={list(obj['grid'])}")
    for k in ("MW", "TPSA", "RotB"):
        g = obj["grid"][k]
        print(f"  {k}: p10={g[10]} p50={g[50]} p90={g[90]}")


if __name__ == "__main__":
    main()
