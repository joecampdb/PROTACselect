"""
Typeset the PROTACselect preprint PDF (reportlab), two-column with a full-width
masthead, embedding the engine-generated figures. Preprint aesthetic.

Run:  <venv>/python report/build_report.py  ->  report/PROTACselect_preprint.pdf
"""
import os

from PIL import Image
from reportlab.lib.colors import HexColor
from reportlab.lib.enums import TA_JUSTIFY, TA_CENTER, TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle, getSampleStyleSheet
from reportlab.lib.units import cm, mm
from reportlab.platypus import (BaseDocTemplate, Frame, PageTemplate, Paragraph,
                                Spacer, Image as RLImage, FrameBreak, NextPageTemplate)

HERE = os.path.dirname(os.path.abspath(__file__))
FIGS = os.path.join(HERE, "figs")
OUT = os.path.join(HERE, "PROTACselect_preprint.pdf")

INK = HexColor("#0f172a"); ACC = HexColor("#4f46e5"); MUTED = HexColor("#64748b")
LINE = HexColor("#d7dbe3"); FAINT = HexColor("#8b93a3")

PW, PHt = A4
LM = RM = 1.55 * cm
TM = 1.5 * cm
BM = 1.5 * cm
GUT = 0.7 * cm
USABLE = PW - LM - RM
COLW = (USABLE - GUT) / 2
MAST_H = 8.9 * cm            # full-width masthead height on page 1

# ---------------------------------------------------------------- styles ----
ss = getSampleStyleSheet()


def S(name, **kw):
    return ParagraphStyle(name, parent=ss["Normal"], **kw)


st_title = S("t", fontName="Helvetica-Bold", fontSize=19, leading=22, textColor=INK,
             alignment=TA_LEFT, spaceAfter=6)
st_auth = S("a", fontName="Helvetica", fontSize=10.5, leading=13, textColor=INK, spaceAfter=2)
st_affil = S("af", fontName="Helvetica-Oblique", fontSize=8.5, leading=11, textColor=MUTED, spaceAfter=1)
st_corr = S("c", fontName="Helvetica", fontSize=8, leading=10, textColor=MUTED, spaceAfter=2)
st_kicker = S("k", fontName="Helvetica-Bold", fontSize=7.5, leading=9, textColor=ACC, spaceAfter=6)
st_abhead = S("abh", fontName="Helvetica-Bold", fontSize=9, leading=11, textColor=INK, spaceBefore=4, spaceAfter=2)
st_abs = S("abs", fontName="Times-Roman", fontSize=8.7, leading=11.6, textColor=INK, alignment=TA_JUSTIFY)
st_kw = S("kw", fontName="Times-Italic", fontSize=8.2, leading=11, textColor=MUTED, spaceBefore=4)
st_h1 = S("h1", fontName="Helvetica-Bold", fontSize=10, leading=12, textColor=INK, spaceBefore=9, spaceAfter=3)
st_body = S("b", fontName="Times-Roman", fontSize=9, leading=11.7, textColor=INK, alignment=TA_JUSTIFY, spaceAfter=5)
st_cap = S("cap", fontName="Helvetica", fontSize=7.5, leading=9.5, textColor=MUTED, spaceBefore=3, spaceAfter=8)
st_ref = S("ref", fontName="Times-Roman", fontSize=7.6, leading=9.4, textColor=INK, spaceAfter=2, leftIndent=10, firstLineIndent=-10)


def figure(fname, colwidth, caption_html):
    path = os.path.join(FIGS, fname)
    iw, ih = Image.open(path).size
    w = colwidth
    h = w * ih / iw
    return [RLImage(path, width=w, height=h), Paragraph(caption_html, st_cap)]


# ---------------------------------------------------------------- content ---
def masthead():
    fw = USABLE
    return [
        Paragraph("PREPRINT · SOFTWARE NOTE · NOT PEER-REVIEWED · 2026-07-08", st_kicker),
        Paragraph("PROTACselect: calibrated physicochemical triage and combinatorial "
                  "exploration of PROTAC design space", st_title),
        Paragraph("J. C.<super>1</super>", st_auth),
        Paragraph("<super>1</super>Independent research · Correspondence: joecamxtc@gmail.com", st_affil),
        Spacer(1, 5),
        Paragraph("ABSTRACT", st_abhead),
        Paragraph(
            "Proteolysis-targeting chimeras (PROTACs) are heterobifunctional degraders assembled from a "
            "target-binding warhead, an E3-ligase ligand, and a connecting linker. The design space is "
            "combinatorial and its productive subset is small, motivating fast triage. We present "
            "<b>PROTACselect</b>, an interactive suite that assembles any warhead·linker·E3 combination, "
            "renders it in 3D, and scores ten physicochemical descriptors against problem-domain objective "
            "functions (weight vectors plus hard feasibility filters). Score thresholds are <i>calibrated</i> "
            "to the physicochemical distribution of 10,728 published degraders from PROTAC-DB, and each metric "
            "is reported with its percentile in that reference set. An Explorer mode sweeps hundreds of "
            "combinations (parametric PEG and alkyl linkers) and ranks them, with a Pareto view of likeness "
            "versus molecular weight. We emphasise, prominently, what the tool does <i>not</i> do: it models "
            "the ligand, not the ternary complex, and its index reflects consistency with known-degrader "
            "chemistry — <i>not</i> validated degradation potency or therapeutic efficacy. PROTACselect is "
            "therefore best understood as an honest developability-triage and design-exploration front-end "
            "for the chemistry half of degrader design.", st_abs),
        Paragraph("<b>Keywords:</b> PROTAC · targeted protein degradation · induced proximity · "
                  "cheminformatics · beyond-Rule-of-5 · design-space exploration", st_kw),
    ]


def body():
    F = []
    F.append(Paragraph("1&nbsp;&nbsp;Introduction", st_h1))
    F.append(Paragraph(
        "PROTACs hijack the ubiquitin–proteasome system: rather than inhibiting a protein, a bifunctional "
        "molecule recruits an E3 ubiquitin ligase to a protein of interest, forming a ternary complex that "
        "tags the target for degradation [1,2]. Because activity is event-driven and catalytic, degraders can "
        "reach targets long considered undruggable [3]. But a PROTAC is defined by three interchangeable "
        "parts, and with even a modest parts list the number of assemblable molecules is large while the "
        "fraction that forms a productive ternary complex is small. Practitioners therefore need rapid triage "
        "that separates the trivially undevelopable from the plausible, before committing to synthesis or to "
        "expensive structure-based modelling.", st_body))
    F.append(Paragraph(
        "PROTACselect addresses the tractable, honest part of this problem — the physicochemistry of the "
        "assembled chimera — while being explicit about the boundary it does not cross.", st_body))

    F += figure("fig1_assembly.png", COLW,
        "<b>Figure 1.</b> Component-to-assembly. A warhead, linker and E3 ligand — each carrying tagged "
        "exit-vector atoms (*:1, *:2) — are fused with RDKit <i>molzip</i> into one degrader. Example: "
        "JQ1 · PEG(n=3) · pomalidomide (CRBN).")

    F.append(Paragraph("2&nbsp;&nbsp;Methods", st_h1))
    F.append(Paragraph(
        "<b>Assembly.</b> Warheads, E3 ligands and linkers are stored as SMILES with tagged exit vectors "
        "(mapped dummy atoms). Assembly joins matching labels via RDKit <i>molzip</i> (Fig. 1). Linkers are "
        "either fixed motifs or generated parametrically — PEG(n=1–12) and alkyl(C2–12).", st_body))
    F.append(Paragraph(
        "<b>Tier-1 descriptors.</b> Ten RDKit descriptors are computed on the assembled molecule: molecular "
        "weight, cLogP, topological polar surface area (TPSA), H-bond donors/acceptors, rotatable bonds, "
        "aromatic rings, fraction sp<super>3</super>, ring count and synthetic accessibility. Each is mapped "
        "through a trapezoidal desirability function to a colour band and combined into a weighted 0–100 index.", st_body))
    F.append(Paragraph(
        "<b>Domains as objectives.</b> A problem domain is an objective function: a weight vector over the ten "
        "descriptors plus optional hard filters. A CNS domain caps TPSA and H-bond donors; a BCL-xL / apoptosis "
        "domain <i>requires</i> a VHL handle, encoding the platelet-sparing rationale of clinical degraders such "
        "as DT2216 [8]. Selecting a domain reweights the index and gates feasibility.", st_body))
    F.append(Paragraph(
        "<b>Calibration.</b> PROTACs occupy beyond-Rule-of-5 space [7], where hand-set drug-likeness cut-offs "
        "mislead. We therefore calibrate every threshold to the empirical descriptor distribution of 10,728 "
        "assembled degraders from PROTAC-DB [4] (Fig. 2): lower-is-safer metrics take green/red edges at the "
        "75th/95th percentiles; the index thus measures consistency with known-degrader space. Each reported "
        "value carries its percentile in this reference set.", st_body))

    F += figure("fig2_space.png", COLW,
        "<b>Figure 2.</b> Physicochemical space of 10,728 PROTAC-DB degraders (dashed line: median). These "
        "distributions define PROTACselect's calibrated bands and per-metric percentile context.")

    F.append(Paragraph("3&nbsp;&nbsp;Results", st_h1))
    F.append(Paragraph(
        "<b>Hand-set thresholds were systematically too strict.</b> Calibration revealed that intuition-based "
        "cut-offs would flag much of real PROTAC space as poor (Fig. 3). The starkest case is TPSA: a hand-set "
        "red edge of 180&nbsp;Å<super>2</super> sits below the median (≈199&nbsp;Å<super>2</super>) of published "
        "degraders. Rotatable bonds and molecular weight show the same bias. Recalibration moves the edges to "
        "where degraders actually live.", st_body))

    F += figure("fig3_calibration.png", COLW,
        "<b>Figure 3.</b> Hand-set “red” edges (dashed) versus data-calibrated 75th/95th-percentile edges "
        "(green/red) on four descriptors. For TPSA and rotatable bonds the hand-set edge lies left of the "
        "bulk — it would reject the majority of real PROTACs.")

    F += figure("fig4_peg.png", COLW,
        "<b>Figure 4.</b> Tier-1 index versus PEG linker length (JQ1·pomalidomide). Developability declines "
        "monotonically with length — a caution: the productive-geometry optimum for <i>degradation</i> is a "
        "Tier-2 property this score does not capture.")

    F.append(Paragraph(
        "<b>Exploration and trade-offs.</b> The Explorer sweeps every linker×E3 combination for a chosen "
        "warhead and domain (310 combinations, &lt;1&nbsp;s) and ranks them. Plotting likeness against molecular "
        "weight exposes the trade-off directly and its Pareto frontier (Fig. 5). Domain hard filters behave as "
        "intended: for a BCL-xL warhead under the apoptosis domain, exactly the VHL combinations pass and all "
        "others are gated out.", st_body))

    F += figure("fig5_pareto.png", COLW,
        "<b>Figure 5.</b> Explorer sweep for enzalutamide under a hormone-driven oral domain: 310 candidates by "
        "domain-weighted likeness versus molecular weight, coloured by verdict (grey = domain-blocked). Dashed: "
        "Pareto frontier.")

    F.append(Paragraph("4&nbsp;&nbsp;Limitations", st_h1))
    F.append(Paragraph(
        "These are first-order, not caveats in fine print. (i) <b>The ligand is modelled, not the biology.</b> "
        "There is no ternary complex, cooperativity, or lysine/ubiquitination geometry — the very features that "
        "determine whether a PROTAC degrades. (ii) <b>The index is typicality, not efficacy.</b> It measures how "
        "consistent a molecule is with known-degrader chemistry; a high-scoring molecule may be biologically "
        "inert. (iii) <b>Nothing is validated against a measured outcome</b> (DC<sub>50</sub>, D<sub>max</sub>, "
        "permeability, oral bioavailability). (iv) Static descriptors such as TPSA miss the chameleonic behaviour "
        "that governs beyond-Rule-of-5 permeability [9]. (v) Advisory tips optimise developability and can, if "
        "followed naively, harm degradation (e.g. shortening a linker). (vi) Tier-2 — degradation potency and "
        "ternary geometry — is out of scope here.", st_body))

    F.append(Paragraph("5&nbsp;&nbsp;Availability", st_h1))
    F.append(Paragraph(
        "PROTACselect is a Python/RDKit engine with a Flask + 3Dmol.js interface. The reference dataset derives "
        "from PROTAC-DB [4], whose terms restrict redistribution; calibration files are regenerated locally. "
        "Roadmap: stratified calibration against oral/clinical degraders; a Tier-2 degradation model reported "
        "with uncertainty; and protein/ternary-complex context.", st_body))

    F.append(Paragraph("References", st_h1))
    refs = [
        "Sakamoto K.M. <i>et al.</i> Protacs: chimeric molecules that target proteins to the SCF complex. "
        "<i>PNAS</i> 98, 8554 (2001).",
        "Winter G.E. <i>et al.</i> Phthalimide conjugation as a strategy for in vivo target protein degradation. "
        "<i>Science</i> 348, 1376 (2015).",
        "Békés M., Langley D.R., Crews C.M. PROTAC targeted protein degraders. <i>Nat. Rev. Drug Discov.</i> "
        "21, 181 (2022).",
        "Weng G. <i>et al.</i> PROTAC-DB: an online database of PROTACs. <i>Nucleic Acids Res.</i> 49, D1381 (2021).",
        "Landrum G. RDKit: open-source cheminformatics. https://www.rdkit.org.",
        "Rego N., Koes D. 3Dmol.js: molecular visualization with WebGL. <i>Bioinformatics</i> 31, 1322 (2015).",
        "DeGoey D.A. <i>et al.</i> Beyond the Rule of 5. <i>J. Med. Chem.</i> 61, 2636 (2018).",
        "Khan S. <i>et al.</i> A selective BCL-X<sub>L</sub> PROTAC degrader (DT2216). <i>Nat. Med.</i> 25, 1938 (2019).",
        "Ermondi G., Vallaro M., Caron G. Degraders early developability assessment. <i>Drug Discov. Today</i> 25, 1585 (2020).",
    ]
    for i, r in enumerate(refs, 1):
        F.append(Paragraph(f"[{i}]&nbsp;&nbsp;{r}", st_ref))
    return F


# ---------------------------------------------------------------- canvas ----
def decorate(canvas, doc):
    canvas.saveState()
    canvas.setFillColor(ACC)
    canvas.rect(0, PHt - 0.32 * cm, PW, 0.32 * cm, stroke=0, fill=1)   # top accent band
    canvas.setFillColor(FAINT)
    canvas.setFont("Helvetica", 7)
    canvas.drawString(LM, 0.85 * cm, "PROTACselect · software preprint · not peer-reviewed")
    canvas.drawRightString(PW - RM, 0.85 * cm, f"{doc.page}")
    canvas.setStrokeColor(LINE)
    canvas.setLineWidth(0.5)
    canvas.line(LM, 1.15 * cm, PW - RM, 1.15 * cm)
    canvas.restoreState()


def build():
    doc = BaseDocTemplate(OUT, pagesize=A4, leftMargin=LM, rightMargin=RM,
                          topMargin=TM, bottomMargin=BM, title="PROTACselect preprint",
                          author="J. C.")
    colH = PHt - TM - BM
    mast = Frame(LM, PHt - TM - MAST_H, USABLE, MAST_H, id="mast",
                 leftPadding=0, rightPadding=0, topPadding=0, bottomPadding=6)
    colH1 = MAST_H  # space consumed above columns on page 1
    p1c1 = Frame(LM, BM, COLW, colH - MAST_H - 0.2 * cm, id="p1c1",
                 leftPadding=0, rightPadding=6, topPadding=0, bottomPadding=0)
    p1c2 = Frame(LM + COLW + GUT, BM, COLW, colH - MAST_H - 0.2 * cm, id="p1c2",
                 leftPadding=6, rightPadding=0, topPadding=0, bottomPadding=0)
    c1 = Frame(LM, BM, COLW, colH, id="c1", leftPadding=0, rightPadding=6, topPadding=0, bottomPadding=0)
    c2 = Frame(LM + COLW + GUT, BM, COLW, colH, id="c2", leftPadding=6, rightPadding=0, topPadding=0, bottomPadding=0)

    doc.addPageTemplates([
        PageTemplate(id="first", frames=[mast, p1c1, p1c2], onPage=decorate),
        PageTemplate(id="rest", frames=[c1, c2], onPage=decorate),
    ])

    story = [NextPageTemplate("rest")]   # page 2+ = two columns, no masthead
    story += masthead()
    story.append(FrameBreak())           # jump masthead -> column 1 (page 1)
    story += body()
    doc.build(story)
    print("wrote", OUT)


if __name__ == "__main__":
    build()
