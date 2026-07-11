# PROTACselect Engine

Assembly + Tier-1 scoring core. Reversible-only, 10×10×10 catalog.

## Setup
```
python -m venv .venv
python -m pip install rdkit
```

## Run the demo
```
# from project root
PYTHONIOENCODING=utf-8 python -m engine.demo
```

## Modules
| File | Role |
|------|------|
| `catalog.py`  | 10 warheads + 10 E3 ligands + 10 linkers, each with a **tagged exit vector** (`[*:1]` warhead side, `[*:2]` E3 side). Base SMILES = PubChem-verified (see `../catalog_smiles.md`). |
| `assembly.py` | `assemble(warhead, linker, e3)` fuses the three fragments at matching map labels via RDKit `molzip`. `validate_catalog()` parses every entry. |
| `metrics.py`  | `score_mol()` = 10 Tier-1 descriptors → per-metric traffic light + desirability + 0–100 index + verdict. `linker_geometry()` = length/flex/rigidity sub-panel. |
| `domains.py`  | 5 problem spaces as **weight vectors + hard filters**. `apply_domain()` re-scores and gates. |
| `demo.py`     | validates catalog, smoke-tests all 1000 combos, scores 5 literature-inspired degraders, and shows the domain lens. |

## API sketch
```python
from engine.assembly import assemble
from engine.metrics import score_mol
from engine.domains import apply_domain
from engine.catalog import warhead_smiles, linker_smiles, e3_smiles, e3_class

mol = assemble(warhead_smiles("JQ1"), linker_smiles("PEG3"), e3_smiles("VH032"))
sc  = score_mol(mol)                      # Tier-1 developability
dom = apply_domain("PS4-Apoptosis", sc, e3_class=e3_class("VH032"))
```

## Status / caveats
- All 30 catalog SMILES parse; all 1000 combinations assemble.
- Tier-1 only (instant, physicochemical). DC50/Dmax/cooperativity/ternary-pose = Tier-2 (not built).
- Exit vectors flagged `Medium`/`Low` in `catalog.py` (MDM2, IAP, VH298, fluoro-Hyp) need chemist
  verification vs PROTAC-DB before production.
- SA score loads from RDKit Contrib if available; otherwise it is omitted from the index.
