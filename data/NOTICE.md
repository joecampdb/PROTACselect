# Data provenance & terms

## What lives here

This directory is intentionally **empty of raw data** in the public repository.

PROTACselect's calibration and percentile context are derived from **PROTAC-DB**
(Weng *et al.*, *Nucleic Acids Research* 2021; Hou group, Zhejiang University), a curated
database of published PROTACs.

## What IS shipped (and why it's OK)

- `engine/calibration.json` — the green/amber/red threshold **edges** (per-metric percentiles).
- `engine/distribution.json` — a 101-point **percentile grid** per descriptor.

These files contain **only aggregate summary statistics** (numbers) computed over the dataset —
**no molecule structures, SMILES, names, or records**. They do not reproduce PROTAC-DB and let the
calibrated features work out of the box.

## What is NOT shipped

- The raw PROTAC SMILES (`data/protac_smiles.smi`) and any PROTAC-DB download (`.xlsx`/`.csv`).
  PROTAC-DB's terms restrict redistribution of the compiled dataset, so it is **git-ignored** and
  must be obtained by each user directly from PROTAC-DB.

## Regenerating the calibration locally (optional)

1. Download the PROTAC set from PROTAC-DB and extract the whole-degrader SMILES, one per line, to
   `data/protac_smiles.smi`.
2. Rebuild the aggregate-statistic files:

   ```bash
   python -m engine.calibrate    data/protac_smiles.smi --source "PROTAC-DB"
   python -m engine.distribution data/protac_smiles.smi --source "PROTAC-DB"
   ```

The app runs **without** these files too — it falls back to built-in heuristic thresholds
(and clearly labels itself "uncalibrated" in that case).

> If you intend to publish this repository, confirm PROTAC-DB's current terms and satisfy yourself
> that shipping aggregate statistics is acceptable for your use; if in doubt, delete the two JSON
> files and let users regenerate them.
