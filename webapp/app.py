"""
PROTACselect webapp — Flask backend.

Serves the single-page GUI and exposes the engine:
  GET  /api/catalog          -> dropdown data (warheads, E3, linkers, domains)
  POST /api/assemble         -> assembled SMILES + 3D molblock + Tier-1 metrics + domain gate

Run (from project root):
  python webapp/app.py
then open http://127.0.0.1:5000
"""
import os
import sys

# make the sibling `engine` package importable
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from flask import Flask, jsonify, render_template, request
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors

from engine.catalog import (WARHEADS, E3_LIGANDS, LINKERS,
                            warhead_smiles, e3_smiles, e3_class, linker_smiles)
from engine.assembly import assemble
from engine.metrics import score_mol, linker_geometry, CALIBRATED, CALIB_META
from engine.domains import apply_domain, DOMAINS
from engine.advisor import generate_tip
from engine.distribution import percentile_of, meta as dist_meta

app = Flask(__name__)
app.config["TEMPLATES_AUTO_RELOAD"] = True   # reflect template edits without a restart

# metric key -> (display label, unit)
METRIC_META = {
    "MW":      ("Molecular weight", "Da"),
    "cLogP":   ("cLogP", ""),
    "TPSA":    ("TPSA", "Å²"),
    "HBD":     ("H-bond donors", ""),
    "HBA":     ("H-bond acceptors", ""),
    "RotB":    ("Rotatable bonds", ""),
    "ArRings": ("Aromatic rings", ""),
    "Fsp3":    ("Fraction sp³", ""),
    "Rings":   ("Total rings", ""),
    "SA":      ("Synthetic access.", ""),
}


# --- parametric linkers -------------------------------------------------------
# PEG and alkyl are generated across lengths; the other 8 stay as fixed motifs.
# Linker identifiers: a catalog key ("Triazole") OR a spec "PEG:n" / "Alkyl:m".
MOTIF_KEYS = ["Alkyl-PEG", "Amide", "Triazole", "Piperazine",
              "Piperidine", "Alkyne", "Phenyl", "Spiro"]


def _peg_smiles(n):    return "[*:1]CC" + "OCC" * n + "[*:2]"   # n ethylene-glycol units
def _alkyl_smiles(m):  return "[*:1]" + "C" * m + "[*:2]"       # m methylenes


def expanded_linkers():
    items = [{"key": k, "name": LINKERS[k][0], "group": "Motif", "rigidity": LINKERS[k][1]}
             for k in MOTIF_KEYS]
    items += [{"key": f"PEG:{n}", "name": f"PEG (n={n})", "group": "PEG · flexible",
               "rigidity": "flexible"} for n in range(1, 13)]
    items += [{"key": f"Alkyl:{m}", "name": f"C{m} alkyl", "group": "Alkyl · flexible",
               "rigidity": "flexible"} for m in range(2, 13)]
    return items


def resolve_linker(spec):
    """spec -> (display_name, tagged_smiles)."""
    if ":" in spec:
        fam, n = spec.split(":")
        n = int(n)
        if fam == "PEG":
            return f"PEG (n={n})", _peg_smiles(n)
        if fam == "Alkyl":
            return f"C{n} alkyl", _alkyl_smiles(n)
    return LINKERS[spec][0], linker_smiles(spec)


def _metric_rows(sc):
    rows = []
    for name, cell in sc["metrics"].items():
        label, unit = METRIC_META[name]
        rows.append({"key": name, "label": label, "unit": unit,
                     "value": cell["value"], "color": cell["color"],
                     "desirability": cell["desirability"],
                     "percentile": percentile_of(name, cell["value"])})
    return rows


def _embed_3d(mol):
    """Return a 3D MOL block (with Hs) or None on failure."""
    try:
        m = Chem.AddHs(mol)
        p = AllChem.ETKDGv3()
        p.randomSeed = 0xf00d
        if AllChem.EmbedMolecule(m, p) != 0:
            if AllChem.EmbedMolecule(m, AllChem.ETKDGv2()) != 0:
                return None
        try:
            AllChem.MMFFOptimizeMolecule(m, maxIters=400)
        except Exception:
            pass
        return Chem.MolToMolBlock(m)
    except Exception:
        return None


@app.route("/")
def index():
    return render_template("index.html")


@app.route("/api/catalog")
def catalog():
    return jsonify({
        "warheads": [{"key": k, "name": v[0], "target": v[1], "confidence": v[4]}
                     for k, v in WARHEADS.items()],
        "e3": [{"key": k, "name": v[0], "e3": v[1], "e3_class": v[2], "confidence": v[5]}
               for k, v in E3_LIGANDS.items()],
        "linkers": expanded_linkers(),
        "domains": [{"key": k, "label": v["label"], "preferred_e3": v["preferred_e3"]}
                    for k, v in DOMAINS.items()],
        "calibration": {"calibrated": CALIBRATED, "meta": CALIB_META},
    })


@app.route("/api/assemble", methods=["POST"])
def api_assemble():
    d = request.get_json(force=True)
    wh, lk, e3 = d["warhead"], d["linker"], d["e3"]
    domain = d.get("domain")

    _lname, lsmi = resolve_linker(lk)
    try:
        mol = assemble(warhead_smiles(wh), lsmi, e3_smiles(e3))
    except Exception as exc:
        return jsonify({"error": f"assembly failed: {exc}"}), 400

    sc = score_mol(mol)
    dom = apply_domain(domain, sc, e3_class=e3_class(e3)) if domain else None
    tip = generate_tip(sc, dom, wh=wh, lk=lk, e3=e3)
    molblock = _embed_3d(mol)

    return jsonify({
        "smiles": Chem.MolToSmiles(mol),
        "mw": round(Descriptors.MolWt(mol), 1),
        "molblock": molblock,
        "embed_ok": molblock is not None,
        "metrics": _metric_rows(sc),
        "index": sc["index"],
        "verdict": sc["verdict"],
        "reds": sc["reds"],
        "ambers": sc["ambers"],
        "linker": linker_geometry(lsmi),
        "domain": dom,
        "tip": tip,
        "source": "assembled",
    })


@app.route("/api/score_molblock", methods=["POST"])
def api_score_molblock():
    """Score a user-uploaded MOL/SDF structure (no component tags available)."""
    d = request.get_json(force=True)
    mb = d.get("molblock", "")
    domain = d.get("domain")

    # SDF may contain multiple records; take the first block
    if "$$$$" in mb:
        mb = mb.split("$$$$")[0]
    mol = Chem.MolFromMolBlock(mb)
    if mol is None:
        return jsonify({"error": "could not parse the uploaded MOL/SDF"}), 400

    sc = score_mol(mol)
    dom = apply_domain(domain, sc, e3_class=None) if domain else None   # E3 unknown
    tip = generate_tip(sc, dom)                                         # metrics-only tip

    # keep the uploaded 3D coords if present, else embed
    conf3d = mol.GetNumConformers() > 0 and mol.GetConformer().Is3D()
    molblock = mb if conf3d else _embed_3d(mol)

    return jsonify({
        "smiles": Chem.MolToSmiles(mol),
        "mw": round(Descriptors.MolWt(mol), 1),
        "molblock": molblock,
        "embed_ok": molblock is not None,
        "metrics": _metric_rows(sc),
        "index": sc["index"],
        "verdict": sc["verdict"],
        "reds": sc["reds"],
        "ambers": sc["ambers"],
        "linker": None,          # no tagged linker in an arbitrary upload
        "domain": dom,
        "tip": tip,
        "source": "uploaded",
    })


@app.route("/api/sweep", methods=["POST"])
def api_sweep():
    """Score ALL linker x E3 combos for one warhead+domain and rank them.

    No 3D embedding here (that's the slow part) -> the full 100-combo sweep is fast.
    Ranked passing-first, then by the domain-weighted score (raw index if no domain).
    """
    d = request.get_json(force=True)
    wh = d["warhead"]
    domain = d.get("domain")

    results = []
    for it in expanded_linkers():
        lspec, lname, lsmi = it["key"], it["name"], resolve_linker(it["key"])[1]
        for e3 in E3_LIGANDS:
            try:
                mol = assemble(warhead_smiles(wh), lsmi, e3_smiles(e3))
            except Exception:
                continue
            sc = score_mol(mol)
            dom = apply_domain(domain, sc, e3_class=e3_class(e3)) if domain else None
            score = dom["weighted_index"] if dom else sc["index"]
            results.append({
                "linker": lspec, "linker_name": lname,
                "e3": e3, "e3_name": E3_LIGANDS[e3][0], "e3_class": e3_class(e3),
                "index": sc["index"], "score": score,
                "verdict": sc["verdict"], "mw": round(Descriptors.MolWt(mol), 1),
                "reds": sc["reds"], "ambers": sc["ambers"],
                "passes": (dom["passes"] if dom else True),
            })

    results.sort(key=lambda r: (r["passes"], r["score"] if r["score"] is not None else -1),
                 reverse=True)
    return jsonify({"warhead": wh, "domain": domain, "n": len(results), "results": results})


if __name__ == "__main__":
    # host 0.0.0.0 + $PORT so this also works on a hosting platform (Render, etc.);
    # locally it's still reachable at http://127.0.0.1:5000
    app.run(host="0.0.0.0", port=int(os.environ.get("PORT", 5000)), debug=False)
