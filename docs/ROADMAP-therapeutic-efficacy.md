# Roadmap: from developability triage to a therapeutic-efficacy version

PROTACselect (v0/v1) scores the **chemistry** of an assembled degrader: is it drug-like, is it typical
of known PROTACs, does it fit a domain's physicochemical objective. That is deliberately the *easy,
honest* half of the problem. This document describes what a **therapeutic-efficacy** version — call it
"PROTACselect-Bio" — would actually have to model, why each piece is hard, and a staged plan to build
it without repeating v0's original sin of dressing a heuristic up as a prediction.

> **Guiding principle.** Every efficacy claim must be tied to a measured outcome, reported *with*
> calibrated uncertainty, and validated on held-out data with leakage-proof splits. If a quantity
> cannot be validated, it is labelled a heuristic — not a prediction.

---

## 1. Why "efficacy" is a different problem

Developability is a property of the **molecule**. Therapeutic efficacy is a property of a
**molecule acting on a biological system**. Between "this compound exists" and "this compound is an
effective medicine" sits a causal chain, and a degrader can fail at every link:

```
warhead+E3 binding → ternary complex → cooperativity → ubiquitination
   → degradation (DC50/Dmax) → catalytic turnover vs. hook effect
   → cell permeability & intracellular exposure → selectivity (proteome-wide)
   → PK / tissue exposure / E3 tissue expression → in-vivo efficacy & therapeutic index
```

v1 touches only the physicochemical inputs to *permeability* and *developability*. Everything to the
left (the biophysics of proximity) and much to the right (systems pharmacology) is unmodelled.

---

## 2. What must be modelled (the target quantities)

| Link in the chain | Quantity to predict | Why it's hard |
|---|---|---|
| Binary binding | K_d(warhead·POI), K_d(E3-ligand·E3) | Warhead K_d often known; still assay-dependent |
| Ternary complex | Formation + **cooperativity α** | Emergent from two protein surfaces; not in the ligand |
| Complex stability | Residence time / half-life | Governs productive ubiquitination window |
| Ubiquitination | Lysine accessibility & transfer efficiency | Depends on ternary **pose** + E2 geometry |
| Degradation | **DC₅₀, Dmax, k_deg** | The actual readout; cell- and target-dependent |
| Kinetics | Catalytic turnover, **hook effect** | Bell-shaped dose response; concentration-dependent |
| Exposure | Passive permeability, efflux, intracellular conc. | bRo5 chameleonicity; poorly captured by 2D descriptors |
| Selectivity | Proteome-wide off-target degradation | Needs global (proteomics) data, not one target |
| Systems | PK, tissue exposure, **E3 tissue expression** | Links to expression atlases + PK/PBPK |
| Outcome | In-vivo efficacy, **therapeutic index** | On-target tox (e.g. BCL-xL in platelets) |

---

## 3. Staged technical plan

Each stage is independently useful and shippable. Effort/uncertainty are relative.

### Stage A — Stratified developability calibration *(low effort, high value)*
**Predicts:** a developability signal that means *developability*, not *typicality*.
**How:** recalibrate the "green" reference against a **developability-relevant subset** (oral/clinical
degraders, or PROTAC-DB rows with favourable measured permeability), so typical-but-non-oral molecules
(e.g. large IV tool compounds) score *lower*, not at ceiling.
**Data:** PROTAC-DB activity/ADME columns; a small curated oral/clinical set.
**Risk:** subset is small and heterogeneous; report as a *relative* score with the reference N shown.
*This is the immediate next step and needs no new modelling machinery.*

### Stage B — Ternary-complex modelling *(high effort, high uncertainty)*
**Predicts:** a plausible ternary pose and a **cooperativity proxy**.
**How:** protein–protein docking constrained by the two ligands (PRosettaC-style), or
AlphaFold-Multimer / diffusion docking with PROTAC restraints; derive a proxy for α from buried
surface area, shape/electrostatic complementarity, and pose-ensemble consistency.
**Data:** the few dozen experimentally solved ternary structures (PDB) for benchmarking.
**Risk:** ternary prediction is an open research problem; poses are ensembles, not single answers.
Ship as *hypotheses with confidence*, never as truth. This is the piece that finally puts the
**proteins** into the tool.

### Stage C — Data-driven degradation model *(medium effort, medium uncertainty)*
**Predicts:** **DC₅₀ / Dmax** (with uncertainty), conditioned on warhead, E3, linker, cell line.
**How:** a graph/relational model over the assembled molecule + ternary features from Stage B
(DeepPROTAC-style), trained on curated activity data; **ensemble or conformal** uncertainty.
**Data:** PROTAC-DB DC₅₀/Dmax with cell-line and E3 context; degradation proteomics screens.
**Risk (the big one):** **data leakage.** Must use scaffold- *and* target-disjoint splits — random
splits vastly overstate performance because analogue series repeat. Assay heterogeneity (different
cell lines, timepoints, readouts) adds noise. Report per-split metrics and calibration curves.

### Stage D — Ubiquitination & selectivity *(high effort)*
**Predicts:** likelihood of productive ubiquitination (surface-lysine mapping on the Stage-B pose);
proteome-wide off-target degradation risk.
**Data:** global degradation proteomics (e.g. expression-/multiplexed-proteomics screens).
**Risk:** selectivity is inherently proteome-scale; this is a research effort in its own right.

### Stage E — Systems pharmacology & therapeutic index *(high effort)*
**Predicts:** tissue exposure and an on-target-toxicity flag by combining PK/PBPK with **E3 tissue
expression** (e.g. expression atlases). Turns the existing domain hard-filters (like "VHL to spare
platelets") from hard-coded rules into data-derived estimates.
**Data:** expression atlases (GTEx/Human Protein Atlas), PK datasets.

---

## 4. Data strategy

- **Sources:** PROTAC-DB (activity), solved ternary complexes (PDB), degradation proteomics screens,
  expression atlases, and — ideally — a proprietary or consortium set of make-and-test measurements.
- **Curation:** harmonise assay type, cell line, timepoint, E3, and units; keep provenance per record.
- **Licensing:** several sources restrict redistribution (as PROTAC-DB does here). Ship **models and
  aggregate statistics**, not raw third-party records. Keep a clear NOTICE per dataset.
- **The core pitfall — leakage:** PROTAC datasets are full of near-duplicate analogue series. Random
  train/test splits produce impressive, meaningless numbers. Mandate **scaffold- and target-disjoint**
  splits, and prefer **prospective** evaluation where possible.

## 5. Validation strategy (non-negotiable)

1. **Retrospective, leakage-proof:** scaffold/target-split held-out sets; report per-split, not a
   single cherry-picked number.
2. **Calibration, not just accuracy:** predicted-vs-measured calibration plots and interval coverage
   (does the 90% interval contain the truth 90% of the time?).
3. **Uncertainty first-class:** every prediction ships with an interval; low-confidence predictions
   are visibly flagged, not hidden.
4. **Prospective where feasible:** a handful of designed-then-measured compounds beats any amount of
   retrospective curve-fitting.
5. **Honest failure modes:** document where the model is blind (novel E3s, novel targets, out-of-domain
   chemistry) rather than extrapolating silently.

## 6. Architecture changes in the app

- **Two-tier UX, clearly separated.** Keep Tier-1 (instant, physicochemical) as-is; add a **Tier-2
  panel** that is explicitly *asynchronous* (a job queue + cache), never implied to be instant.
- **Confidence everywhere.** Tier-2 numbers render with intervals and a confidence badge; a low-
  confidence DC₅₀ looks different from a high-confidence one.
- **Mode switch.** A visible toggle between *"developability only"* and *"predicted efficacy"* so users
  never confuse the two — the failure mode v1 was careful to avoid.
- **Domains extended.** A domain objective can then include predicted efficacy *with its uncertainty*,
  not just physicochemistry — e.g. maximise predicted Dmax subject to CNS caps and a platelet-safety
  filter derived from Stage E rather than hard-coded.

## 7. Honest expectations

Even a good therapeutic-efficacy version would be a **prioritisation aid, not a replacement for
experiment.** Degradation prediction is unsolved; ternary modelling is uncertain; efficacy data are
scarce and noisy. The right ambition is a tool that **reorders a synthesis queue** and **flags the
implausible** better than physicochemistry alone — while stating its uncertainty out loud. Anything
that promises to tell you a molecule *will* work would be repeating exactly the mistake this project
was built to avoid.

---

## Summary sequencing

| Stage | Deliverable | Effort | Confidence it can be built well |
|---|---|---|---|
| A | Stratified developability score | Low | High |
| B | Ternary pose + cooperativity proxy | High | Medium |
| C | DC₅₀/Dmax model + uncertainty | Medium | Medium (given clean splits) |
| D | Ubiquitination + selectivity | High | Low–Medium |
| E | PK / tissue / therapeutic index | High | Low–Medium |

Start with **A** (it makes the current score honest again), then **B** (it finally adds the proteins),
then **C** (the first genuine efficacy signal) — each shipped with uncertainty and validated on
leakage-proof splits.
