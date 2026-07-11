# Limitations

PROTACselect is an **honest developability-triage and design-exploration tool for the chemistry half
of degrader design.** It is **not** a predictor of degradation or therapeutic efficacy. These
limitations are first-order — they define what the tool is for — not fine print.

## 1. It models the ligand, not the biology

A PROTAC works by forming a **ternary complex** (target protein · PROTAC · E3 ligase) that positions
the target for ubiquitination and proteasomal degradation. PROTACselect models **none** of this:
there is no protein, no ternary complex, no cooperativity (α), no surface-lysine or ubiquitination
geometry. The 3D viewer shows a single ligand in a vacuum. The "larger structures" that decide
whether a degrader works — the proteins — are absent.

## 2. The score is *typicality*, not efficacy

After calibration to PROTAC-DB, the 0–100 index measures **how consistent a molecule is with
known-degrader physicochemical space** — i.e. "does this look like a real PROTAC?". It does **not**
measure "will this degrade the target." A molecule can score highly and be biologically inert.

## 3. Nothing is validated against a measured outcome

No score here has been regressed against DC₅₀, Dmax, cellular degradation, permeability, or oral
bioavailability. The calibration grounds the *thresholds* in real chemistry; it does **not** make the
index an outcome model.

## 4. The descriptors are imperfect in beyond-Rule-of-5 space

PROTACs are large and flexible. Calculated logP is extrapolated well outside its training regime, and
a static TPSA misses the **molecular chameleonicity** (dynamic intramolecular H-bonding) that actually
governs beyond-Rule-of-5 permeability. The descriptors are useful signals, not ground truth.

## 5. A single 3D conformer is not meaningful

The viewer embeds one gas-phase conformer of a molecule with 15–30 rotatable bonds and an enormous
conformational ensemble. It is illustrative, not the bioactive or ternary-bound conformation.

## 6. Most of the combinatorial grid is non-functional — and the tool scores it anyway

The overwhelming majority of assemblable combinations would not degrade anything. PROTACselect scores
each one confidently on developability; it **cannot** distinguish "developable and active" from
"developable and inert."

## 7. Advisory tips optimise developability and can point the wrong way

The Tip advisor improves the (developability) score. Some of that advice — e.g. *shorten the linker to
cut molecular weight* — can **destroy** degradation, because linker length vs. activity is
non-monotonic and the shortest linker often breaks the productive ternary geometry. Tips are
developability-only and should be read as such.

## 8. Some catalog entries carry known uncertainty

Several exit-vector placements and SMILES in `engine/catalog.py` are flagged Medium/Low confidence
(notably MDM2/IAP handles, VH298, and a fluorinated-hydroxyproline VHL ligand). Assembly forms bonds
without a chemical-sanity check, so some assembled molecules are synthetically implausible. Verify
against primary literature before serious use.

## 9. Tier-2 is out of scope

Degradation potency (DC₅₀, Dmax), cooperativity, ternary-complex geometry, selectivity, and PK/tissue
exposure are **not** modelled. Bridging that gap is a research programme, sketched in
[ROADMAP-therapeutic-efficacy.md](ROADMAP-therapeutic-efficacy.md).

---

**In one sentence:** use PROTACselect to quickly triage and explore degrader *chemistry*; do **not**
use it to conclude that a molecule will degrade a target or be therapeutically effective.
