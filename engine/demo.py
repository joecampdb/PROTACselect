"""
PROTACselect engine demo.

Run:  python -m engine.demo     (from project root)
"""
import sys
try:
    sys.stdout.reconfigure(encoding="utf-8")   # emoji-safe on Windows consoles
except Exception:
    pass

from rdkit import Chem
from rdkit.Chem import Descriptors

from .catalog import (WARHEADS, E3_LIGANDS, LINKERS,
                      warhead_smiles, e3_smiles, e3_class, linker_smiles)
from .assembly import assemble, validate_catalog
from .metrics import score_mol, linker_geometry
from .domains import apply_domain, DOMAINS

DOT = {"green": "🟢", "amber": "🟡", "red": "🔴", "n/a": "⚪"}


def build_and_score(wh, lk, e3):
    mol = assemble(warhead_smiles(wh), linker_smiles(lk), e3_smiles(e3))
    return mol, score_mol(mol)


def print_score(title, wh, lk, e3, mol, sc):
    print(f"\n=== {title} ===")
    print(f"  {wh} + {lk} + {e3}   (MW {Descriptors.MolWt(mol):.0f})")
    print(f"  SMILES: {Chem.MolToSmiles(mol)}")
    cells = []
    for name, c in sc["metrics"].items():
        v = c["value"]
        vs = f"{v:.1f}" if isinstance(v, float) else str(v)
        cells.append(f"{DOT[c['color']]}{name}={vs}")
    print("  " + "  ".join(cells))
    lg = linker_geometry(linker_smiles(lk))
    print(f"  linker: len={lg['length_atoms']}({DOT[lg['length_band']]}) "
          f"rotB={lg['rot_bonds']}({DOT[lg['rot_band']]}) "
          f"Fsp3={lg['fsp3']}({DOT[lg['rigid_band']]})")
    print(f"  --> index={sc['index']}  reds={sc['reds']} ambers={sc['ambers']}  "
          f"VERDICT: {sc['verdict']}"
          + ("" if sc["sa_available"] else "   [SA score unavailable]"))


def main():
    print("### Catalog validation ###")
    fails = validate_catalog(WARHEADS, E3_LIGANDS, LINKERS)
    if fails:
        for kind, key, smi in fails:
            print(f"  FAIL {kind}:{key}  {smi}")
    else:
        n = len(WARHEADS) + len(E3_LIGANDS) + len(LINKERS)
        print(f"  all {n} catalog SMILES parse OK "
              f"({len(WARHEADS)} warheads, {len(E3_LIGANDS)} E3, {len(LINKERS)} linkers)")

    # exhaustive assembly smoke test: every warhead x linker x E3
    print("\n### Full-grid assembly smoke test (10x10x10) ###")
    ok = err = 0
    errors = []
    for wh in WARHEADS:
        for lk in LINKERS:
            for e3 in E3_LIGANDS:
                try:
                    assemble(warhead_smiles(wh), linker_smiles(lk), e3_smiles(e3))
                    ok += 1
                except Exception as e:
                    err += 1
                    errors.append((wh, lk, e3, str(e)))
    print(f"  assembled {ok}/{ok + err} combinations")
    for wh, lk, e3, msg in errors[:10]:
        print(f"  ERR {wh}+{lk}+{e3}: {msg}")

    # representative literature-inspired degraders
    examples = [
        ("MZ1-like  (BET/VHL)",      "JQ1",         "PEG3",       "VH032"),
        ("dBET-like (BET/CRBN)",     "JQ1",         "Alkyl-C6",   "Pomalidomide"),
        ("ARV-471-like (ERa/CRBN)",  "4-OHT",       "Amide",      "Pomalidomide"),
        ("ARV-110-like (AR/CRBN)",   "Enzalutamide","Piperidine", "Pomalidomide"),
        ("DT2216-like (BCLxL/VHL)",  "A-1155463",   "PEG3",       "VH032"),
    ]
    for title, wh, lk, e3 in examples:
        mol, sc = build_and_score(wh, lk, e3)
        print_score(title, wh, lk, e3, mol, sc)

    # domain application: score DT2216-like under every domain (shows hard filters)
    print("\n### Domain lens: DT2216-like across all problem spaces ###")
    mol, sc = build_and_score("A-1155463", "PEG3", "VH032")
    for dk in DOMAINS:
        res = apply_domain(dk, sc, e3_class=e3_class("VH032"))
        gate = "PASS" if res["passes"] else "BLOCKED"
        fr = "; ".join(f"{k}:{'ok' if v else 'X'}" for k, v in res["filter_results"].items()) or "none"
        print(f"  [{gate:7}] {dk:22} idx={res['weighted_index']:<5} filters: {fr}")


if __name__ == "__main__":
    main()
