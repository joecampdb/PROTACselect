"""
PROTACselect catalog: warheads, E3 ligands, linkers.

Each ligand SMILES carries a TAGGED EXIT VECTOR as a mapped dummy atom:
  - warhead exit  -> [*:1]   (bonds to the linker's [*:1])
  - E3 exit       -> [*:2]   (bonds to the linker's [*:2])
  - linker        -> both [*:1] and [*:2]

The assembly engine joins matching map numbers via RDKit molzip:
    warhead[*:1] -- [*:1]linker[*:2] -- [*:2]E3

NOTE ON PROVENANCE: base SMILES were sourced verbatim from PubChem (see catalog_smiles.md
for CIDs). The [*:n] placement encodes the exit vector described there. Exit vectors marked
confidence != "High" (esp. MDM2/IAP, VH298, fluoro-Hyp) are best-effort and flagged for a
chemist to verify against PROTAC-DB before production use.
"""

# ---------------------------------------------------------------- warheads ---
# key: (name, target, tagged_smiles, exit_vector_desc, confidence)
WARHEADS = {
    "JQ1":         ("(+)-JQ1", "BET/BRD4",
                    "CC1=C(SC2=C1C(=N[C@H](C3=NN=C(N32)C)CC(=O)[*:1])C4=CC=C(C=C4)Cl)C",
                    "tert-butyl ester -> COOH, amide-couple linker at carbonyl", "High"),
    "Dasatinib":   ("Dasatinib", "BCR-ABL/SRC",
                    "CC1=C(C(=CC=C1)Cl)NC(=O)C2=CN=C(S2)NC3=CC(=NC(=N3)C)N4CCN(CC4)[*:1]",
                    "piperazine N (hydroxyethyl removed)", "High"),
    "Palbociclib": ("Palbociclib", "CDK4/6",
                    "CC1=C(C(=O)N(C2=NC(=NC=C12)NC3=NC=C(C=C3)N4CCN([*:1])CC4)C5CCCC5)C(=O)C",
                    "terminal piperazine N-H", "High"),
    "Gefitinib":   ("Gefitinib", "EGFR",
                    "COC1=C(C=C2C(=C1)N=CN=C2NC3=CC(=C(C=C3)F)Cl)OCCC[*:1]",
                    "morpholino-propoxy tail (morpholine replaced)", "Medium"),
    "Lapatinib":   ("Lapatinib", "HER2/EGFR",
                    "CS(=O)(=O)CCN([*:1])CC1=CC=C(O1)C2=CC3=C(C=C2)N=CN=C3NC4=CC(=C(C=C4)OCC5=CC(=CC=C5)F)Cl",
                    "(methylsulfonyl)ethyl-amino tail", "Medium"),
    "Ceritinib":   ("Ceritinib", "ALK",
                    "CC1=CC(=C(C=C1C2CCN([*:1])CC2)OC(C)C)NC3=NC=C(C(=N3)NC4=CC=CC=C4S(=O)(=O)C(C)C)Cl",
                    "piperidine N-H (MS4078/TL13-112)", "High"),
    "Vemurafenib": ("Vemurafenib", "BRAF V600E",
                    "CCCS(=O)(=O)N([*:1])C1=C(C(=C(C=C1)F)C(=O)C2=CNC3=C2C=C(C=N3)C4=CC=C(C=C4)Cl)F",
                    "propyl-sulfonamide N", "Medium-Low"),
    "Enzalutamide":("Enzalutamide", "Androgen Receptor",
                    "CC1(C(=O)N(C(=S)N1C2=CC(=C(C=C2)C(=O)N[*:1])F)C3=CC(=C(C=C3)C#N)C(F)(F)F)C",
                    "N-methyl carboxamide (ARCC-4; note ARV-110 differs)", "High"),
    "4-OHT":       ("4-Hydroxytamoxifen (Z)", "ERalpha",
                    "CC/C(=C(\\C1=CC=C(C=C1)O)/C2=CC=C(C=C2)OCCN([*:1])C)/C3=CC=CC=C3",
                    "dimethylaminoethoxy tail (note ARV-471 uses lasofoxifene)", "Medium"),
    "A-1155463":   ("A-1155463", "BCL-xL",
                    "[*:1]CC#CC1=CC(=C(C=C1)OCCCC2=C(N=C(S2)N3CCC4=C(C3)C(=CC=C4)C(=O)NC5=NC6=CC=CC=C6S5)C(=O)O)F",
                    "dimethylamino-butynyl terminus (note DT2216 uses navitoclax)", "Medium"),
}

# --------------------------------------------------------------- E3 ligands ---
# key: (name, e3, e3_class, tagged_smiles, exit_vector_desc, confidence)
E3_LIGANDS = {
    "Pomalidomide":("Pomalidomide", "CRBN", "CRBN",
                    "C1CC(=O)NC(=O)C1N2C(=O)C3=C(C2=O)C(=CC=C3)N[*:2]",
                    "aromatic 4-amino (canonical CRBN vector)", "High"),
    "Lenalidomide":("Lenalidomide", "CRBN", "CRBN",
                    "C1CC(=O)NC(=O)C1N2CC3=C(C2=O)C=CC=C3N[*:2]",
                    "aromatic 4-amino", "High"),
    "Thalidomide": ("Thalidomide (4-O)", "CRBN", "CRBN",
                    "C1CC(=O)NC(=O)C1N2C(=O)C3=CC=CC(=C3C2=O)O[*:2]",
                    "phthalimide C4 via 4-hydroxy ether", "High"),
    "VH032":       ("VH032/(S,R,S)-AHPC", "VHL", "VHL",
                    "CC1=C(SC=N1)C2=CC=C(C=C2)CNC(=O)[C@@H]3C[C@H](CN3C(=O)[C@H](C(C)(C)C)NC(=O)[*:2])O",
                    "terminal N-acetyl cap (MZ1/ACBI)", "High"),
    "VH298":       ("VH298", "VHL", "VHL",
                    "CC1=C(SC=N1)C2=CC=C(C=C2[*:2])CNC(=O)[C@@H]3C[C@H](CN3C(=O)[C@H](C(C)(C)C)NC(=O)C4(CC4)C#N)O",
                    "benzyl phenyl (meta) -- placeholder, verify", "Medium"),
    "FluoroHyp":   ("Fluoro-Hyp VHL (Ciulli 14a)", "VHL", "VHL",
                    "CC1=C(SC=N1)C2=CC=C(C=C2)CNC(=O)[C@@H]3[C@H](F)[C@H](CN3C(=O)[C@H](C(C)(C)C)NC(=O)[*:2])O",
                    "N-acetyl cap; SMILES DERIVED - verify vs PDB 6GFY", "Low"),
    "Nutlin-3a":   ("Nutlin-3a", "MDM2", "MDM2",
                    "CC(C)OC1=C(C=CC(=C1)OC)C2=N[C@H]([C@H](N2C(=O)N3CCN([*:2])C(=O)C3)C4=CC=C(C=C4)Cl)C5=CC=C(C=C5)Cl",
                    "piperazinone NH (A1874-like) -- verify", "Medium"),
    "Idasanutlin": ("Idasanutlin/RG7388", "MDM2", "MDM2",
                    "CC(C)(C)C[C@H]1[C@]([C@H]([C@@H](N1)C(=O)NC2=C(C=C(C=C2)C(=O)[*:2])OC)C3=C(C(=CC=C3)Cl)F)(C#N)C4=C(C=C(C=C4)Cl)F",
                    "terminal COOH -- verify", "Medium"),
    "AMG232":      ("AMG 232/navtemadlin", "MDM2", "MDM2",
                    "CC(C)[C@@H](CS(=O)(=O)C(C)C)N1[C@@H]([C@H](C[C@](C1=O)(C)CC(=O)[*:2])C2=CC(=CC=C2)Cl)C3=CC=C(C=C3)Cl",
                    "terminal acetic-acid COOH", "Medium-High"),
    "LCL161":      ("LCL161", "IAP (cIAP/XIAP)", "IAP",
                    "C[C@@H](C(=O)N[C@@H](C1CCCCC1)C(=O)N2CCC[C@H]2C3=NC(=CS3)C(=O)C4=CC=C(C=C4)[*:2])NC",
                    "4-fluorobenzoyl para (F replaced) -- verify", "Medium-Low"),
}

# ------------------------------------------------------------------ linkers ---
# key: (name, class, rigidity, tagged_smiles)
LINKERS = {
    "PEG3":       ("PEG (n=3)", "flexible", "[*:1]CCOCCOCCOCC[*:2]"),
    "Alkyl-C6":   ("Hexyl", "flexible", "[*:1]CCCCCC[*:2]"),
    "Alkyl-PEG":  ("Alkyl-PEG hybrid", "flexible", "[*:1]CCCOCCOCC[*:2]"),
    "Amide":      ("Amide-containing", "flexible-semi", "[*:1]CCC(=O)NCCC[*:2]"),
    "Triazole":   ("1,2,3-Triazole (click)", "semi-rigid", "[*:1]CCn1cc(CC[*:2])nn1"),
    "Piperazine": ("Piperazine", "semi-rigid", "[*:1]N1CCN(CC1)[*:2]"),
    "Piperidine": ("Piperidine (1,4)", "semi-rigid", "[*:1]N1CCC(CC1)[*:2]"),
    "Alkyne":     ("Alkyne rod", "rigid", "[*:1]CC#CC[*:2]"),
    "Phenyl":     ("p-Phenylene", "rigid", "[*:1]c1ccc(cc1)[*:2]"),
    "Spiro":      ("Spiro[3.3]heptane", "rigid-3D", "[*:1]C1CC2(C1)CC([*:2])C2"),
}


def warhead_smiles(key):  return WARHEADS[key][2]
def e3_smiles(key):       return E3_LIGANDS[key][3]
def e3_class(key):        return E3_LIGANDS[key][2]
def linker_smiles(key):   return LINKERS[key][2]
