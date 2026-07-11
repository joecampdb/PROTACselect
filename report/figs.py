"""
Generate publication figures for the PROTACselect preprint — straight from the
engine + PROTAC-DB data (not GUI screenshots), so every figure is reproducible.

Run:  <venv>/python report/figs.py
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

from engine.metrics import METRICS, score_mol
from engine.catalog import (warhead_smiles, e3_smiles, linker_smiles,
                            E3_LIGANDS, e3_class)
from engine.assembly import assemble
from engine.domains import apply_domain

HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.join(HERE, "figs")
os.makedirs(FIGS, exist_ok=True)
DATA = os.path.join(os.path.dirname(HERE), "data", "protac_smiles.smi")

INK = "#0f172a"; ACC = "#4f46e5"; ACC2 = "#7c3aed"
GREEN = "#16a34a"; AMBER = "#e08a00"; RED = "#dc2626"
MUTED = "#64748b"; GRID = "#e6e9ee"; FILL = "#c7cdf2"
VC = {"Typical": GREEN, "Borderline": AMBER, "Atypical": RED}

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 9,
    "axes.edgecolor": "#cbd2dc", "axes.linewidth": 0.8,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.7, "axes.axisbelow": True,
    "xtick.color": MUTED, "ytick.color": MUTED, "axes.labelcolor": INK, "text.color": INK,
    "figure.dpi": 220, "savefig.dpi": 220,
})

UNITS = {"MW": "Da", "cLogP": "", "TPSA": "Å²", "HBD": "", "HBA": "",
         "RotB": "", "ArRings": "", "Fsp3": "", "Rings": "", "SA": ""}
ORDER = ["MW", "cLogP", "TPSA", "HBD", "HBA", "RotB", "ArRings", "Fsp3", "Rings", "SA"]


def peg(n):   return "[*:1]CC" + "OCC" * n + "[*:2]"
def alkyl(m): return "[*:1]" + "C" * m + "[*:2]"


def load_descriptors():
    smis = [s for s in open(DATA, encoding="utf-8", errors="ignore").read().split() if s]
    cols = {k: [] for k in METRICS}
    n = 0
    for s in smis:
        m = Chem.MolFromSmiles(s)
        if m is None:
            continue
        n += 1
        for k, (fn, _b, _w, _h) in METRICS.items():
            try:
                v = fn(m)
                if v is not None:
                    cols[k].append(float(v))
            except Exception:
                pass
    return {k: np.array(v) for k, v in cols.items()}, n


# --- Fig 1: component -> assembly schematic ----------------------------------
def fig_assembly():
    def img(smi, size):
        return Draw.MolToImage(Chem.MolFromSmiles(smi), size=size)
    wh = img(warhead_smiles("JQ1"), (300, 210))
    lk = img(linker_smiles("PEG3"), (300, 150))
    e3 = img(e3_smiles("Pomalidomide"), (300, 210))
    prod = Draw.MolToImage(
        assemble(warhead_smiles("JQ1"), linker_smiles("PEG3"), e3_smiles("Pomalidomide")),
        size=(680, 240))

    fig = plt.figure(figsize=(7.1, 4.5))
    gs = fig.add_gridspec(2, 3, height_ratios=[1, 1.15], hspace=0.28, wspace=0.05)
    for ax, im, title, col in [
        (fig.add_subplot(gs[0, 0]), wh, "Warhead · JQ1 (BET)", "#2563eb"),
        (fig.add_subplot(gs[0, 1]), lk, "Linker · PEG (n=3)", "#0d9488"),
        (fig.add_subplot(gs[0, 2]), e3, "E3 ligand · Pomalidomide (CRBN)", "#9333ea"),
    ]:
        ax.imshow(im); ax.axis("off")
        ax.set_title(title, fontsize=8.5, fontweight="bold", color=col, pad=4)
    for x in (0.345, 0.66):
        fig.text(x, 0.72, "+", ha="center", va="center", fontsize=15, color=MUTED, fontweight="bold")
    axp = fig.add_subplot(gs[1, :])
    axp.imshow(prod); axp.axis("off")
    axp.set_title("↓  assembled PROTAC  ·  molzip join at tagged exit vectors",
                  fontsize=8.5, fontweight="bold", color=INK, pad=6)
    fig.savefig(os.path.join(FIGS, "fig1_assembly.png"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


# --- Fig 2: PROTAC-DB physicochemical space ----------------------------------
def fig_space(desc, n):
    fig, axes = plt.subplots(2, 5, figsize=(7.1, 3.2))
    for ax, k in zip(axes.flat, ORDER):
        a = desc[k]
        lo, hi = np.percentile(a, 0.5), np.percentile(a, 99.0)   # clip extreme outliers
        if hi <= lo:
            hi = lo + 1
        ax.hist(a, bins=36, range=(lo, hi), color=ACC, alpha=0.85, edgecolor="white", linewidth=0.25)
        ax.axvline(np.median(a), color=INK, lw=1.0, ls="--")
        ax.set_xlim(lo, hi)
        u = f" ({UNITS[k]})" if UNITS[k] else ""
        ax.set_title(f"{k}{u}", fontsize=8, fontweight="bold")
        ax.set_yticks([]); ax.tick_params(labelsize=6.5)
    fig.tight_layout(pad=0.6)
    fig.savefig(os.path.join(FIGS, "fig2_space.png"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


# --- Fig 3: invented thresholds vs data-calibrated ---------------------------
INVENT = {"MW": 1050, "TPSA": 180, "RotB": 16, "HBD": 6}   # original heuristic "red" edge


def fig_calibration(desc):
    keys = ["MW", "TPSA", "RotB", "HBD"]
    fig, axes = plt.subplots(1, 4, figsize=(7.1, 2.5))
    for ax, k in zip(axes, keys):
        a = desc[k]
        p75, p95 = np.percentile(a, 75), np.percentile(a, 95)
        hi = max(np.percentile(a, 98), INVENT[k], p95) * 1.12       # keep both edges + bulk visible
        ax.hist(a, bins=34, range=(0, hi), color=FILL, edgecolor="white", linewidth=0.25)
        ax.axvline(INVENT[k], color=MUTED, lw=1.5, ls="--")
        ax.axvline(p75, color=GREEN, lw=1.4)
        ax.axvline(p95, color=RED, lw=1.6)
        ax.set_xlim(0, hi)
        u = f" ({UNITS[k]})" if UNITS[k] else ""
        ax.set_title(f"{k}{u}", fontsize=8.5, fontweight="bold")
        ax.set_yticks([]); ax.tick_params(labelsize=6.5)
    handles = [plt.Line2D([], [], color=MUTED, ls="--", lw=1.5, label="hand-set “red” edge"),
               plt.Line2D([], [], color=GREEN, lw=1.4, label="calibrated P75 (green)"),
               plt.Line2D([], [], color=RED, lw=1.6, label="calibrated P95 (red)")]
    fig.legend(handles=handles, ncol=3, fontsize=7, frameon=False,
               loc="lower center", bbox_to_anchor=(0.5, -0.06))
    fig.tight_layout(pad=0.6)
    fig.savefig(os.path.join(FIGS, "fig3_calibration.png"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


# --- Fig 4: PEG length vs PROTAC-likeness (developability, not degradation) --
def fig_peg():
    wh, e3 = warhead_smiles("JQ1"), e3_smiles("Pomalidomide")
    ns = list(range(1, 13))
    scores, verdicts = [], []
    for nn in ns:
        sc = score_mol(assemble(wh, peg(nn), e3))
        scores.append(sc["index"]); verdicts.append(sc["verdict"])
    fig, ax = plt.subplots(figsize=(3.45, 3.0))
    ax.plot(ns, scores, "-", color=ACC, lw=1.8, zorder=1)
    ax.scatter(ns, scores, c=[VC[v] for v in verdicts], s=42, zorder=2,
               edgecolors="white", linewidths=0.8)
    ax.set_xlabel("PEG linker length (n glycol units)")
    ax.set_ylabel("PROTAC-likeness index")
    ax.set_ylim(0, 105); ax.set_xticks(ns[::2])
    ax.annotate("shorter ⇒ higher Tier-1 score\n(≠ higher degradation)",
                xy=(9, scores[8]), xytext=(5.4, 32), fontsize=6.8, color=MUTED,
                arrowprops=dict(arrowstyle="->", color=MUTED, lw=0.8))
    fig.tight_layout(pad=0.5)
    fig.savefig(os.path.join(FIGS, "fig4_peg.png"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


# --- Fig 5: Pareto sweep (likeness vs MW) ------------------------------------
def fig_pareto():
    wh = warhead_smiles("Enzalutamide"); domain = "PS2-HormoneOral"
    specs = [linker_smiles(k) for k in
             ["Alkyl-PEG", "Amide", "Triazole", "Piperazine", "Piperidine", "Alkyne", "Phenyl", "Spiro"]]
    specs += [peg(n) for n in range(1, 13)] + [alkyl(m) for m in range(2, 13)]
    pts = []
    for ls in specs:
        for e3k in E3_LIGANDS:
            mol = assemble(wh, ls, e3_smiles(e3k))
            sc = score_mol(mol)
            dom = apply_domain(domain, sc, e3_class=e3_class(e3k))
            pts.append((Descriptors.MolWt(mol), dom["weighted_index"],
                        dom["passes"], sc["verdict"]))
    fig, ax = plt.subplots(figsize=(3.45, 3.0))
    for mw, s, passes, verd in pts:
        ax.scatter(mw, s, s=26, c=(VC[verd] if passes else MUTED),
                   alpha=(0.9 if passes else 0.22),
                   edgecolors="white", linewidths=0.5, zorder=2)
    # Pareto frontier: maximize likeness, minimize MW
    ok = sorted([(mw, s) for mw, s, p, v in pts], key=lambda t: (t[0], -t[1]))
    front, best = [], -1
    for mw, s in ok:
        if s > best:
            front.append((mw, s)); best = s
    fx, fy = [m for m, _ in front], [s for _, s in front]
    ax.plot(fx, fy, color=ACC2, lw=2.0, ls="--", zorder=4, label="Pareto frontier")
    ax.scatter(fx, fy, s=20, color=ACC2, edgecolors="white", linewidths=0.6, zorder=5)
    ax.set_xlabel("Molecular weight (Da)")
    ax.set_ylabel("Domain-weighted likeness")
    ax.set_ylim(0, 105)
    ax.legend(fontsize=7, frameon=False, loc="lower left")
    fig.tight_layout(pad=0.5)
    fig.savefig(os.path.join(FIGS, "fig5_pareto.png"), bbox_inches="tight", facecolor="white")
    plt.close(fig)


if __name__ == "__main__":
    print("loading descriptors from PROTAC-DB set…")
    desc, n = load_descriptors()
    print(f"  {n} molecules")
    fig_assembly();     print("  fig1 assembly")
    fig_space(desc, n); print("  fig2 space")
    fig_calibration(desc); print("  fig3 calibration")
    fig_peg();          print("  fig4 peg")
    fig_pareto();       print("  fig5 pareto")
    print("done ->", FIGS)
