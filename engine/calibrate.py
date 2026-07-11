"""
Tier-1 threshold calibration.

Turns the hand-set green/amber/red bands in metrics.py into DATA-DRIVEN bands
derived from the descriptor distribution of real published PROTACs, and writes
engine/calibration.json (which metrics.load_calibration() then applies).

Per-metric directionality is preserved:
  lower-is-safer (MW,TPSA,HBD,HBA,RotB,ArRings,SA): green <= P75, amber P75–P95, red > P95
  higher-is-better (Fsp3):                          green >= P25, amber P05–P25, red < P05
  two-sided (cLogP, Rings):                         green P25–P75, amber to P05/P95, red beyond

IMPORTANT (see log.md Limitations): this calibrates against *typical PROTAC space*
— "are you consistent with known degraders?" — NOT validated efficacy or oral-%F.
It grounds the thresholds in data; it does not make the index an outcome model.

CLI:
  <venv>/python -m engine.calibrate data/protac_smiles.smi --source "PROTAC-DB"
"""
import argparse
import json
import math
import os

import numpy as np
from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")
from .metrics import METRICS  # reuse the exact compute fns + directionality flags

INF = math.inf


def _ser(x):
    if x == INF:
        return "inf"
    if x == -INF:
        return "-inf"
    return round(float(x), 3)


def read_smiles(path, smiles_col=None):
    """Read SMILES from .smi/.txt (first token/line) or .csv/.tsv (a smiles column)."""
    out = []
    lower = path.lower()
    delim = "," if lower.endswith(".csv") else ("\t" if lower.endswith(".tsv") else None)
    with open(path, encoding="utf-8", errors="ignore") as fh:
        first = True
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if delim:
                parts = line.split(delim)
                if first:
                    first = False
                    header = [p.strip().lower() for p in parts]
                    if smiles_col is None:
                        for i, h in enumerate(header):
                            if "smiles" in h:
                                smiles_col = i
                                break
                    if smiles_col is not None:   # skip header row
                        continue
                idx = smiles_col if smiles_col is not None else 0
                if idx < len(parts):
                    out.append(parts[idx].strip().strip('"'))
            else:
                out.append(line.split()[0])
                first = False
    return [s for s in out if s and s.lower() != "smiles"]


def _row(mol):
    r = {}
    for k, (fn, _b, _w, _h) in METRICS.items():
        try:
            r[k] = fn(mol)
        except Exception:
            r[k] = None
    return r


def calibrate(smiles, source, min_valid=50):
    rows, n_read, n_ok = [], 0, 0
    for smi in smiles:
        n_read += 1
        m = Chem.MolFromSmiles(smi)
        if m is None:
            continue
        rows.append(_row(m))
        n_ok += 1

    bands, pct = {}, {}
    for k, (fn, old, w, hib) in METRICS.items():
        vals = np.array([r[k] for r in rows if r.get(k) is not None], dtype=float)
        if vals.size < min_valid:
            bands[k] = [_ser(v) for v in old]      # keep heuristic if too little data
            continue
        P = lambda q: float(np.percentile(vals, q))
        if hib is True:               # higher is better
            b = (P(5), P(25), INF, INF)
        elif hib is False:            # lower is safer
            b = (-INF, -INF, P(75), P(95))
        else:                         # two-sided
            b = (P(5), P(25), P(75), P(95))
        bands[k] = [_ser(v) for v in b]
        pct[k] = {str(q): round(P(q), 2) for q in (5, 25, 50, 75, 95)}

    return {"meta": {"source": source, "n_read": n_read, "n_valid": n_ok},
            "bands": bands, "percentiles": pct}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("smiles_file")
    ap.add_argument("--source", default="unknown")
    ap.add_argument("--smiles-col", type=int, default=None)
    ap.add_argument("--out", default=os.path.join(os.path.dirname(__file__), "calibration.json"))
    a = ap.parse_args()

    smis = read_smiles(a.smiles_file, a.smiles_col)
    print(f"read {len(smis)} SMILES from {a.smiles_file}")
    res = calibrate(smis, a.source)
    with open(a.out, "w") as fh:
        json.dump(res, fh, indent=2)
    print(f"n_valid={res['meta']['n_valid']}  -> wrote {a.out}\n")
    for k in METRICS:
        print(f"  {k:8} bands={res['bands'][k]}  pct={res['percentiles'].get(k)}")


if __name__ == "__main__":
    main()
