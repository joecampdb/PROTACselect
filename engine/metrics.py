"""
Tier-1 metric panel (see tier1_metrics.md).

Instant, RDKit-native, PROTAC (beyond-Rule-of-5) calibrated descriptors + a
weighted 0-100 developability index and a traffic-light verdict.

Each metric declares a desirability trapezoid (rl, gl, gh, rh):
    d = 1.0        for gl <= x <= gh          (green)
    d ramps 1->0   across [rl,gl] and [gh,rh] (amber)
    d = 0.0        for x <= rl or x >= rh      (red)
Use +/-inf for one-sided metrics.
"""
import math
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Crippen

INF = math.inf

# --- SA score is an RDKit Contrib module; load it if present -----------------
try:
    import os, sys
    from rdkit.Chem import RDConfig
    sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
    import sascorer  # type: ignore
    _HAVE_SA = True
except Exception:
    _HAVE_SA = False

# metric -> (compute_fn, (rl, gl, gh, rh), weight, higher_is_better)
METRICS = {
    "MW":       (lambda m: Descriptors.MolWt(m),               (-INF, -INF, 900, 1050), 0.14, False),
    "cLogP":    (lambda m: Crippen.MolLogP(m),                 (1, 2, 5, 7),            0.12, None),
    "TPSA":     (lambda m: Descriptors.TPSA(m),                (-INF, -INF, 140, 180),  0.16, False),
    "HBD":      (lambda m: rdMolDescriptors.CalcNumHBD(m),     (-INF, -INF, 3, 6),      0.18, False),
    "HBA":      (lambda m: rdMolDescriptors.CalcNumHBA(m),     (-INF, -INF, 10, 15),    0.08, False),
    "RotB":     (lambda m: rdMolDescriptors.CalcNumRotatableBonds(m), (-INF, -INF, 10, 16), 0.10, False),
    "ArRings":  (lambda m: rdMolDescriptors.CalcNumAromaticRings(m),  (-INF, -INF, 4, 6),   0.06, False),
    "Fsp3":     (lambda m: rdMolDescriptors.CalcFractionCSP3(m),(0.2, 0.3, INF, INF),    0.08, True),
    "Rings":    (lambda m: rdMolDescriptors.CalcNumRings(m),    (1.5, 3, 6, 8),          0.03, None),
    "SA":       (lambda m: sascorer.calculateScore(m) if _HAVE_SA else None,
                                                               (-INF, -INF, 4, 6),      0.05, False),
}


# --- optional data-driven calibration ---------------------------------------
# If engine/calibration.json exists (produced by calibrate.py from a real PROTAC
# dataset), its percentile-derived bands REPLACE the hand-set heuristics above.
# Weights are NOT changed here (they remain editorial — see log.md Limitations).
import json as _json

_CALIB_PATH = os.path.join(os.path.dirname(__file__), "calibration.json")
CALIBRATED = False
CALIB_META = None


def _parse_band(v):
    if v == "inf":
        return INF
    if v == "-inf":
        return -INF
    return float(v)


def load_calibration(path=_CALIB_PATH):
    """Overlay percentile-derived bands from calibration.json onto METRICS."""
    global CALIBRATED, CALIB_META
    if not os.path.exists(path):
        return False
    with open(path) as fh:
        data = _json.load(fh)
    for key, band in data.get("bands", {}).items():
        if key in METRICS:
            fn, _old, w, hib = METRICS[key]
            METRICS[key] = (fn, tuple(_parse_band(x) for x in band), w, hib)
    CALIBRATED = True
    CALIB_META = data.get("meta")
    return True


load_calibration()


def _desirability(x, rl, gl, gh, rh):
    if x is None:
        return None
    if x < gl:
        if rl == -INF:
            return 1.0
        return 0.0 if x <= rl else (x - rl) / (gl - rl)
    if x > gh:
        if rh == INF:
            return 1.0
        return 0.0 if x >= rh else (rh - x) / (rh - gh)
    return 1.0


def _color(x, rl, gl, gh, rh):
    if x is None:
        return "n/a"
    if (rl != -INF and x <= rl) or (rh != INF and x >= rh):
        return "red"
    if x < gl or x > gh:
        return "amber"
    return "green"


def score_mol(mol, weights=None):
    """Return dict: per-metric {value,color,desirability} + composite + verdict."""
    w = dict((k, METRICS[k][2]) for k in METRICS)
    if weights:
        w.update(weights)

    per = {}
    wsum = 0.0
    dsum = 0.0
    reds = ambers = 0
    for name, (fn, bands, _, _) in METRICS.items():
        val = fn(mol)
        des = _desirability(val, *bands)
        col = _color(val, *bands)
        per[name] = {"value": val, "color": col, "desirability": des}
        if col == "red":
            reds += 1
        elif col == "amber":
            ambers += 1
        if des is not None:                      # skip missing SA from the index
            dsum += w[name] * des
            wsum += w[name]

    index = round(100 * dsum / wsum, 1) if wsum else None
    # verdicts describe consistency with typical PROTAC space (post-calibration),
    # NOT validated developability/efficacy — see log.md Limitations.
    if reds >= 2:
        verdict = "Atypical"
    elif reds == 1 or ambers >= 3:
        verdict = "Borderline"
    else:
        verdict = "Typical"

    return {"metrics": per, "index": index, "reds": reds,
            "ambers": ambers, "verdict": verdict, "sa_available": _HAVE_SA}


# ------------------------------------------------------- linker sub-panel ----
def linker_geometry(linker_smiles):
    """Length / rotatable bonds / rigidity of a bis-[*] linker fragment."""
    m = Chem.MolFromSmiles(linker_smiles)
    if m is None:
        return None
    dummies = [a.GetIdx() for a in m.GetAtoms() if a.GetAtomicNum() == 0]
    if len(dummies) != 2:
        return {"error": f"expected 2 attachment points, found {len(dummies)}"}
    # heavy-atom neighbours of the two dummies = the real attachment atoms
    a1 = m.GetAtomWithIdx(dummies[0]).GetNeighbors()[0].GetIdx()
    a2 = m.GetAtomWithIdx(dummies[1]).GetNeighbors()[0].GetIdx()
    path = Chem.GetShortestPath(m, a1, a2)
    length = len(path)                      # atoms spanning the two anchors, inclusive
    rotb = rdMolDescriptors.CalcNumRotatableBonds(m)
    fsp3 = rdMolDescriptors.CalcFractionCSP3(m)
    has_ring = rdMolDescriptors.CalcNumRings(m) > 0

    def band(v, lo_ok, hi_ok, lo_edge, hi_edge):
        if v < lo_edge or v > hi_edge:
            return "red"
        if v < lo_ok or v > hi_ok:
            return "amber"
        return "green"

    return {
        "length_atoms": length,   "length_band": band(length, 6, 16, 4, 20),
        "rot_bonds": rotb,        "rot_band": band(rotb, 2, 8, 1, 12),
        "fsp3": round(fsp3, 2),   "rigid_band": "green" if (fsp3 >= 0.4 or has_ring) else
                                                ("amber" if fsp3 >= 0.2 else "red"),
    }
