"""
Advisor: for a low-scoring design, produce an actionable Tip that accounts for
the chosen problem domain AND the 3-component structure (warhead / linker / E3).

Triggered when the effective developability score (domain-weighted if a domain
is selected, else the raw Tier-1 index) is < 50.0.
"""
from .catalog import WARHEADS, E3_LIGANDS, LINKERS

# linker character (by catalog key)
RIGID   = {"Alkyne", "Phenyl", "Spiro"}
LONG    = {"PEG3", "Alkyl-PEG", "Alkyl-C6"}
POLAR   = {"PEG3", "Alkyl-PEG", "Amide"}      # add TPSA / H-bond donors

TIP_FLOOR = 50.0


def _num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def generate_tip(score_result, domain_result=None, wh=None, lk=None, e3=None):
    """Return a 'Tip: ...' string, or None if the design scores >= 50."""
    idx = score_result.get("index")
    score = (domain_result.get("weighted_index")
             if domain_result and domain_result.get("weighted_index") is not None
             else idx)
    if score is None or score >= TIP_FLOOR:
        return None

    metrics = score_result["metrics"]
    reds = {k for k, c in metrics.items() if c["color"] == "red"}
    val = {k: _num(c["value"]) for k, c in metrics.items()}

    e3cls = E3_LIGANDS[e3][2] if e3 else None
    e3name = E3_LIGANDS[e3][0] if e3 else None

    clauses = []

    # 1) size ----------------------------------------------------------------
    if "MW" in reds:
        bits = []
        if e3cls in ("VHL", "MDM2", "IAP"):
            bits.append("pair it with a smaller CRBN handle (Pomalidomide/Lenalidomide) to shed ~100–200 Da")
        if lk in LONG:
            bits.append("shorten the linker (try the Alkyne rod or p-Phenylene)")
        joiner = " and ".join(bits) if bits else "trim the linker and pick a smaller E3 handle"
        clauses.append((3, f"MW ({val['MW']:.0f}) is over the bRo5 ceiling — {joiner}."))

    # 2) flexibility ---------------------------------------------------------
    if "RotB" in reds and lk not in RIGID:
        clauses.append((3, f"{int(val['RotB'])} rotatable bonds is very flexible — a rigid linker "
                            "(p-Phenylene, Alkyne, or Spiro[3.3]heptane) pre-organizes the ternary "
                            "complex and improves permeability."))

    # 3) polarity too high ---------------------------------------------------
    if ("TPSA" in reds or "HBD" in reds):
        extra = ""
        if lk in POLAR:
            extra = f" The {LINKERS[lk][0]} linker adds polar surface/donors — "
        clauses.append((2, f"Polar surface / H-bond donors are high.{extra} an alkyl or aryl linker "
                            "lowers TPSA and donor count."))

    # 4) lipophilicity -------------------------------------------------------
    if "cLogP" in reds and val["cLogP"] is not None:
        if val["cLogP"] > 7:
            clauses.append((2, f"cLogP ({val['cLogP']:.1f}) is too greasy — a PEG linker adds polarity to rebalance."))
        elif val["cLogP"] < 1:
            clauses.append((2, f"cLogP ({val['cLogP']:.1f}) is too polar to permeate — an alkyl linker restores permeability."))

    # 5) aromatic overload ---------------------------------------------------
    if "ArRings" in reds:
        clauses.append((1, f"{int(val['ArRings'])} aromatic rings — pick a saturated linker "
                            "(Spiro[3.3]heptane or hexyl) over Phenyl/Triazole."))

    # 6) domain fit ----------------------------------------------------------
    if domain_result and not domain_result.get("passes", True):
        failed = [k for k, v in domain_result.get("filter_results", {}).items() if v is False]
        if any("E3 class" in f for f in failed) and e3cls and e3cls != "VHL":
            clauses.append((4, f"This domain needs a VHL handle to spare platelets — swap the E3 to "
                                f"VH032; {e3name} won't gate through."))
        if any(f.split()[0] in ("TPSA", "HBD", "MW", "cLogP") for f in failed):
            clauses.append((4, "These CNS caps are strict — if brain penetration isn't essential, re-frame "
                                "under a peripheral (PS3a) or hematologic (PS1) domain that tolerates bRo5 size."))
    elif not domain_result:
        clauses.append((1, "Under an oral/hormone domain this would score even lower; a hematologic "
                           "(PS1-Heme) framing tolerates this size better."))

    if not clauses:
        clauses.append((1, "Several descriptors sit in the red band — try the shortest, most rigid linker "
                           "and the smallest E3 handle your domain allows."))

    # highest-severity first, keep it to 2 clauses
    clauses.sort(key=lambda c: c[0], reverse=True)
    return "Tip: " + " ".join(text for _, text in clauses[:2])
