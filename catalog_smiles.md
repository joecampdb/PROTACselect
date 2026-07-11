# PROTACselect — Catalog SMILES + Exit Vectors (sourced 2026-07-06)

> SMILES pulled live from **PubChem PUG-REST** this session (copied verbatim). PubChem CID is the
> verified anchor. `[*]` = tagged exit-vector atom placeholder (where the linker bonds). Exit-vector
> SMILES below are **annotations to still be encoded** — the engine must place `[*]` at the named atom.
> Confidence reflects how standardized the exit vector is across published degraders.

## Legend
- **CID** = PubChem CID (verified this session).
- **Exit vector** = solvent-exposed atom/position where linker attaches; binding-critical groups are noted as OFF-LIMITS.
- ⚠️ = needs verification before hard-coding SMILES/stereo.

---

## Warheads (target-binding ligands)

| # | Name | Target | CID | SMILES (PubChem, verbatim) | Exit vector | Conf |
|---|------|--------|-----|-----------------------------|-------------|------|
| 1 | (+)-JQ1 | BET/BRD4 | 46907787 | `CC1=C(SC2=C1C(=N[C@H](C3=NN=C(N32)C)CC(=O)OC(C)(C)C)C4=CC=C(C=C4)Cl)C` | **tert-butyl ester → hydrolyze to COOH, amide-couple** (dBET1, MZ1, ARV-825). t-Bu ester points to solvent. | High |
| 2 | Dasatinib | BCR-ABL/SRC | 3062316 | `CC1=C(C(=CC=C1)Cl)NC(=O)C2=CN=C(S2)NC3=CC(=NC(=N3)C)N4CCN(CC4)CCO` | **piperazine N** (remove hydroxyethyl, attach directly) or terminal **–OH**. | High |
| 3 | Palbociclib | CDK4/6 | 5330286 | `CC1=C(C(=O)N(C2=NC(=NC=C12)NC3=NC=C(C=C3)N4CCNCC4)C5CCCC5)C(=O)C` | **terminal piperazine N–H**. | High |
| 4 | Gefitinib | EGFR | 123631 | `COC1=C(C=C2C(=C1)N=CN=C2NC3=CC(=C(C=C3)F)Cl)OCCCN4CCOCC4` | **morpholino-propoxy tail** (replace morpholine / derivatize alkoxy). Less standardized. | Med |
| 5 | Lapatinib | HER2/EGFR | 208908 | `CS(=O)(=O)CCNCC1=CC=C(O1)C2=CC3=C(C=C2)N=CN=C3NC4=CC(=C(C=C4)OCC5=CC(=CC=C5)F)Cl` | **(methylsulfonyl)ethyl-amino tail** (secondary amine / sulfonyl terminus). | Med |
| 6 | Ceritinib (LDK378) | ALK | 57379345 | `CC1=CC(=C(C=C1C2CCNCC2)OC(C)C)NC3=NC=C(C(=N3)NC4=CC=CC=C4S(=O)(=O)C(C)C)Cl` | **piperidine N–H** (MS4078, TL13-12/112). | High |
| 7 | Vemurafenib | BRAF V600E | 42611257 | `CCCS(=O)(=O)NC1=C(C(=C(C=C1)F)C(=O)C2=CNC3=C2C=C(C=N3)C4=CC=C(C=C4)Cl)F` | **n-propyl-sulfonamide** region. Few degraders → less canonical. | Med-Low |
| 8 | Enzalutamide | Androgen Receptor | 15951529 | `CC1(C(=O)N(C(=S)N1C2=CC(=C(C=C2)C(=O)NC)F)C3=CC(=C(C=C3)C#N)C(F)(F)F)C` | **N-methyl carboxamide** (replace N-methyl) — ARCC-4 vector. ⚠️ ARV-110 uses a DIFFERENT proprietary AR ligand. | High |
| 9 | 4-Hydroxytamoxifen (Z) | ERα | 449459 | `CC/C(=C(\C1=CC=C(C=C1)O)/C2=CC=C(C=C2)OCCN(C)C)/C3=CC=CC=C3` | **dimethylaminoethoxy tail** (replace amine); phenol = alt. ⚠️ ARV-471 uses **lasofoxifene** (CID 216416), exit = pyrrolidinyl-ethoxy tail. | Med |
| 10 | A-1155463 | BCL-xL | 59447577 | `CN(C)CC#CC1=CC(=C(C=C1)OCCCC2=C(N=C(S2)N3CCC4=C(C3)C(=CC=C4)C(=O)NC5=NC6=CC=CC=C6S5)C(=O)O)F` | **dimethylamino-butynyl (alkyne) terminus**. ⚠️ Clinical DT2216 uses **navitoclax** (CID 24978538), linker at aryl-**piperazine**, recruits VHL. | High (navitoclax vector) |

## E3 ligase ligands (recruiters)

| # | Name | E3 | CID | SMILES (PubChem, verbatim) | Exit vector | Conf |
|---|------|----|-----|-----------------------------|-------------|------|
| 1 | Pomalidomide | CRBN | 134780 | `C1CC(=O)NC(=O)C1N2C(=O)C3=C(C2=O)C(=CC=C3)N` | **aromatic 4-amino (C4-NH₂)** → amide/amine. Canonical CRBN vector. Alt: glutarimide N–H (rare). | High |
| 2 | Lenalidomide | CRBN | 216326 | `C1CC(=O)NC(=O)C1N2CC3=C(C2=O)C=CC=C3N` | **aromatic 4-amino (C4-NH₂)**. Alt: glutarimide N–H. | High |
| 3 | Thalidomide | CRBN | 5426 | `C1CC(=O)NC(=O)C1N2C(=O)C3=CC=CC=C3C2=O` | **phthalimide C4** via 4-F/4-OH/4-NH₂ analogs (ether/SNAr/amine). Parent lacks C4 substituent. | High |
| 4 | VH032 / (S,R,S)-AHPC | VHL | 77232228 | `CC1=C(SC=N1)C2=CC=C(C=C2)CNC(=O)[C@@H]3C[C@H](CN3C(=O)[C@H](C(C)(C)C)NC(=O)C)O` | (a) **terminal N-acetyl cap** (MZ1/ACBI classic); (b) **benzyl phenyl ring** (para/meta). Hydroxyproline-OH + thiazole OFF-LIMITS. | High |
| 5 | VH298 | VHL | 122199236 | `CC1=C(SC=N1)C2=CC=C(C=C2)CNC(=O)[C@@H]3C[C@H](CN3C(=O)[C@H](C(C)(C)C)NC(=O)C4(CC4)C#N)O` | Same positions as VH032 (cap region / benzyl phenyl). Mostly used as free VHL inhibitor. | Med |
| 6 | Fluoro-Hyp VHL variant (3F,4OH-Pro; Ciulli 2018 "14a") | VHL | ⚠️ none clean | ⚠️ **UNVERIFIED** (derived): `CC1=C(SC=N1)C2=CC=C(C=C2)CNC(=O)[C@@H]3[C@H](F)[C@H](CN3C(=O)[C@H](C(C)(C)C)NC(=O)C)O` — PDB 6GFY. Verify C3-F stereo. | Same as VH032 (F = metabolic mod, NOT linker point). | Struct High / SMILES Low |
| 7 | Nutlin-3a | MDM2 | 11433190 | `CC(C)OC1=C(C=CC(=C1)OC)C2=N[C@H]([C@H](N2C(=O)N3CCNC(=O)C3)C4=CC=C(C=C4)Cl)C5=CC=C(C=C5)Cl` | **piperazinone / methoxy "eastern" side chain** (A1874). Chlorophenyls + imidazoline OFF-LIMITS. | Med |
| 8 | Idasanutlin (RG7388) | MDM2 | 53358942 | `CC(C)(C)C[C@H]1[C@]([C@H]([C@@H](N1)C(=O)NC2=C(C=C(C=C2)C(=O)O)OC)C3=C(C(=CC=C3)Cl)F)(C#N)C4=C(C=C(C=C4)Cl)F` | **terminal COOH** of the methoxy-carboxy-phenyl arm (amide couple). | Med |
| 9 | AMG 232 / navtemadlin | MDM2 | 58573469 | `CC(C)[C@@H](CS(=O)(=O)C(C)C)N1[C@@H]([C@H](C[C@](C1=O)(C)CC(=O)O)C2=CC(=CC=C2)Cl)C3=CC=C(C=C3)Cl` | **terminal COOH** (acetic-acid arm on piperidinone quaternary C). | Med-High |
| 10 | LCL161 | IAP (cIAP/XIAP) | 24737642 | `C[C@@H](C(=O)N[C@@H](C1CCCCC1)C(=O)N2CCC[C@H]2C3=NC(=CS3)C(=O)C4=CC=C(C=C4)F)NC` | **4-fluorobenzoyl / thiazole "eastern" region**. N-methyl-Ala AVPI N-terminus OFF-LIMITS (BIR-binding). | Med-Low |

---

## Caveats (carry forward)
- **ChEMBL IDs deliberately omitted / unverified** — use CID to cross-link. Only A-1155463 = CHEMBL3342332 was cross-checked.
- **Row 6 (fluoro-Hyp VHL) SMILES is DERIVED, not PubChem-verified** — check against PDB 6GFY (esp. C3-F stereocenter) before hard-coding.
- **Stereochemistry:** JQ1 = (S)/(+); 4-OHT = (Z)/afimoxifene (CID 449459; 5352135 = (E)); navitoclax/lasofoxifene isomeric SMILES retain stereocenters.
- **Exit-vector confidence:** CRBN 4-NH₂ and VHL cap/phenyl = High (crystal-guided, many degraders). MDM2 + IAP = Med/Low (fewer degraders, series-dependent) → confirm per-series in PROTAC-DB (protac-db.org) before building.
- **Clinical-vs-classic warhead flags:** ARV-110 ≠ enzalutamide; ARV-471 = lasofoxifene (not 4-OHT); DT2216 = navitoclax+VHL (not A-1155463+CRBN).
- **Next:** encode `[*]` at each named exit atom; define the 10 linker chemotypes as bis-`[*]` SMILES; write the assembly routine (join warhead-`[*]` + `[*]`-linker-`[*]` + `[*]`-E3).
