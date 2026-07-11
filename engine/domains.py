"""
Problem domains = the objective function. Each domain is:
  - a WEIGHT VECTOR over the Tier-1 metrics (reweights the developability index), and
  - optional HARD FILTERS (metric caps and/or a mandatory E3 class) that gate feasibility.

Derived from warhead -> target -> disease biology (see log.md).
"""

# op helpers for hard filters on metric values
def _le(v): return ("<=", v)
def _ge(v): return (">=", v)
def _between(lo, hi): return ("between", (lo, hi))

DOMAINS = {
    "PS1-Heme": {
        "label": "Hematologic / epigenetic (marrow, no BBB)",
        "warheads": ["JQ1", "Dasatinib"],
        "weights": {},                         # balanced default
        "hard_filters": {},
        "preferred_e3": "CRBN",
    },
    "PS2-HormoneOral": {
        "label": "Hormone-driven solid tumor (chronic oral)",
        "warheads": ["Enzalutamide", "4-OHT", "Palbociclib"],
        "weights": {"HBD": 0.22, "TPSA": 0.20, "MW": 0.18, "cLogP": 0.12,
                    "RotB": 0.10, "HBA": 0.06, "Fsp3": 0.04, "ArRings": 0.04,
                    "SA": 0.02, "Rings": 0.02},
        "hard_filters": {},
        "preferred_e3": "CRBN",
    },
    "PS3a-KinasePeripheral": {
        "label": "Driver-kinase solid tumor (peripheral)",
        "warheads": ["Gefitinib", "Lapatinib", "Ceritinib", "Vemurafenib"],
        "weights": {},
        "hard_filters": {},
        "preferred_e3": "CRBN",
    },
    "PS3b-KinaseCNS": {
        "label": "Driver-kinase solid tumor (CNS-penetrant, brain mets)",
        "warheads": ["Gefitinib", "Ceritinib", "Vemurafenib"],
        "weights": {"TPSA": 0.26, "HBD": 0.24, "MW": 0.16, "cLogP": 0.12,
                    "RotB": 0.08, "ArRings": 0.05, "Fsp3": 0.04, "HBA": 0.02,
                    "SA": 0.02, "Rings": 0.01},
        "hard_filters": {"TPSA": _le(90), "HBD": _le(2), "MW": _le(750),
                         "cLogP": _between(1, 4)},
        "preferred_e3": "CRBN",
    },
    "PS4-Apoptosis": {
        "label": "Apoptosis / BCL-xL (toxicity-gated, spare platelets)",
        "warheads": ["A-1155463"],
        "weights": {},
        "hard_filters": {"__e3_class__": "VHL"},   # mechanistic requirement (DT2216)
        "preferred_e3": "VHL",
    },
}


def _check(op, target, value):
    if value is None:
        return False
    if op == "<=":
        return value <= target
    if op == ">=":
        return value >= target
    if op == "between":
        lo, hi = target
        return lo <= value <= hi
    return False


def apply_domain(domain_key, score_result, e3_class=None):
    """
    Re-score with the domain weight vector and evaluate hard filters.
    Returns {index, weighted_index, passes, filter_results}.
    """
    d = DOMAINS[domain_key]
    per = score_result["metrics"]

    # weighted index under this domain's weights
    from .metrics import METRICS
    w = {k: METRICS[k][2] for k in METRICS}
    w.update(d["weights"])
    dsum = wsum = 0.0
    for name, cell in per.items():
        des = cell["desirability"]
        if des is not None:
            dsum += w[name] * des
            wsum += w[name]
    weighted_index = round(100 * dsum / wsum, 1) if wsum else None

    # hard filters (ok True/False, or None = indeterminate -> does not block)
    filter_results = {}
    passes = True
    for key, spec in d["hard_filters"].items():
        if key == "__e3_class__":
            ok = None if e3_class is None else (e3_class == spec)
            filter_results[f"E3 class == {spec}"] = ok
        else:
            op, target = spec
            ok = _check(op, target, per.get(key, {}).get("value"))
            filter_results[f"{key} {op} {target}"] = ok
        if ok is False:
            passes = False

    return {"domain": domain_key, "label": d["label"],
            "weighted_index": weighted_index, "passes": passes,
            "filter_results": filter_results}
