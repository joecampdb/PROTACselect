"use strict";

const $ = (id) => document.getElementById(id);
const COLORS = { green: "#16a34a", amber: "#e08a00", red: "#dc2626", "n/a": "#9aa4b2" };
const GAUGE_LEN = 2 * Math.PI * 52; // r=52

let CAT = null;
let viewer = null;
let curStyle = "stick";
let lastMolblock = null;

// ---------------------------------------------------------------- boot ------
async function boot() {
  CAT = await (await fetch("/api/catalog")).json();

  fillSelect($("warhead"), CAT.warheads, (o) => `${o.name} — ${o.target}`);
  fillLinkerSelect(CAT.linkers);
  fillSelect($("e3"), CAT.e3, (o) => `${o.name} — ${o.e3}`);
  fillSelect($("domain"), CAT.domains, (o) => o.label, true);

  // show calibration provenance under the gauge
  const cal = CAT.calibration;
  if (cal && cal.calibrated && cal.meta) {
    $("index-sub").textContent = `Calibrated to ${cal.meta.source || "PROTAC-DB"} · n=${(cal.meta.n_valid || 0).toLocaleString()}`;
  } else {
    $("index-sub").textContent = "Uncalibrated heuristic thresholds";
  }

  // sensible defaults (an oral-degrader-shaped combo)
  $("warhead").value = "Enzalutamide";
  $("linker").value = "Piperidine";
  $("e3").value = "Pomalidomide";

  viewer = $3Dmol.createViewer($("viewer"), { backgroundColor: 0x0a0e1c });

  ["warhead", "linker", "e3", "domain"].forEach((id) =>
    $(id).addEventListener("change", update));

  $("style-seg").addEventListener("click", (e) => {
    const b = e.target.closest("button"); if (!b) return;
    document.querySelectorAll("#style-seg button").forEach((x) => x.classList.remove("active"));
    b.classList.add("active"); curStyle = b.dataset.style; applyStyle();
  });
  $("copy").addEventListener("click", () => {
    navigator.clipboard.writeText($("smiles").textContent).then(() => {
      const b = $("copy"); b.textContent = "Copied"; setTimeout(() => (b.textContent = "Copy"), 1200);
    });
  });

  // download current 3D structure as SDF
  $("download").addEventListener("click", () => {
    if (!lastMolblock) return;
    const sdf = lastMolblock.replace(/\s*$/, "") + "\n$$$$\n";
    const a = document.createElement("a");
    a.href = URL.createObjectURL(new Blob([sdf], { type: "chemical/x-mdl-sdfile" }));
    a.download = "protacselect.sdf"; a.click(); URL.revokeObjectURL(a.href);
  });

  // upload a MOL/SDF and score it
  $("upload-btn").addEventListener("click", () => $("upload-file").click());
  $("upload-file").addEventListener("change", uploadMol);

  // explorer drawer
  $("explore-btn").addEventListener("click", openExplore);
  $("drawer-close").addEventListener("click", closeDrawer);
  $("drawer-backdrop").addEventListener("click", closeDrawer);
  document.addEventListener("keydown", (e) => { if (e.key === "Escape") closeDrawer(); });

  update();
}

// ------------------------------------------------------------- explorer -----
function openDrawer() { $("drawer").classList.add("open"); $("drawer-backdrop").classList.add("show"); }
function closeDrawer() { $("drawer").classList.remove("open"); $("drawer-backdrop").classList.remove("show"); }

async function openExplore() {
  const wh = $("warhead").value;
  const domain = $("domain").value || null;
  const whName = metaFor(CAT.warheads, wh).name;
  const domLabel = domain ? metaFor(CAT.domains, domain).label : "no domain";
  $("drawer-title").textContent = "Top candidates";
  $("drawer-sub").textContent = `${whName} · ${domLabel}`;
  $("drawer-list").innerHTML = '<div class="drawer-empty">Ranking 100 combos…</div>';
  openDrawer();
  let data;
  try {
    data = await (await fetch("/api/sweep", {
      method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ warhead: wh, domain }),
    })).json();
  } catch (e) { $("drawer-list").innerHTML = '<div class="drawer-empty">Sweep failed.</div>'; return; }
  renderSweep(data);
}

function loadCombo(linker, e3) {
  $("linker").value = linker;
  $("e3").value = e3;
  closeDrawer();
  update();
}

function renderScatter(results) {
  const host = $("drawer-scatter");
  const pts = results.filter((r) => r.score != null);
  if (!pts.length) { host.innerHTML = ""; return; }
  const W = 448, H = 210, padL = 40, padR = 14, padT = 14, padB = 30;
  const mws = pts.map((p) => p.mw);
  let minX = Math.min(...mws), maxX = Math.max(...mws);
  if (minX === maxX) { minX -= 10; maxX += 10; }
  const X = (v) => padL + ((v - minX) / (maxX - minX)) * (W - padL - padR);
  const Y = (s) => H - padB - (s / 100) * (H - padT - padB);

  // Pareto front: maximize likeness, minimize MW
  const sorted = [...pts].sort((a, b) => a.mw - b.mw || b.score - a.score);
  let best = -1; const front = [];
  sorted.forEach((p) => { if (p.score > best) { front.push(p); best = p.score; } });
  const frontLine = front.map((p) => `${X(p.mw).toFixed(1)},${Y(p.score).toFixed(1)}`).join(" ");

  const vmap = { Typical: "green", Borderline: "amber", Atypical: "red" };
  const grid = [0, 50, 100].map((t) =>
    `<line class="axl" x1="${padL}" y1="${Y(t)}" x2="${W - padR}" y2="${Y(t)}"/>` +
    `<text class="ax" x="${padL - 6}" y="${Y(t) + 3}" text-anchor="end">${t}</text>`).join("");
  const dots = pts.map((p) => {
    const c = p.passes ? (vmap[p.verdict] || "na") : "blocked";
    return `<circle class="sdot ${c}" cx="${X(p.mw).toFixed(1)}" cy="${Y(p.score).toFixed(1)}" r="4"` +
      ` data-linker="${p.linker}" data-e3="${p.e3}">` +
      `<title>${p.e3_name} + ${p.linker_name}\nlikeness ${p.score} · ${p.mw} Da</title></circle>`;
  }).join("");

  host.innerHTML =
    `<svg viewBox="0 0 ${W} ${H}" class="scatter">
       ${grid}
       ${front.length > 1 ? `<polyline class="front" points="${frontLine}"/>` : ""}
       ${dots}
       <text class="axt" x="10" y="${padT + 1}">likeness ↑</text>
       <text class="ax" x="${padL}" y="${H - padB + 15}" text-anchor="middle">${Math.round(minX)}</text>
       <text class="axt" x="${(padL + W - padR) / 2}" y="${H - 5}" text-anchor="middle">MW (Da) →</text>
       <text class="ax" x="${W - padR}" y="${H - padB + 15}" text-anchor="end">${Math.round(maxX)}</text>
     </svg>`;
  host.querySelectorAll(".sdot").forEach((d) =>
    d.addEventListener("click", () => loadCombo(d.dataset.linker, d.dataset.e3)));
}

function renderSweep(data) {
  renderScatter(data.results);
  const vmap = { Typical: "green", Borderline: "amber", Atypical: "red" };
  const host = $("drawer-list");
  host.innerHTML = "";
  data.results.forEach((r, i) => {
    const vc = vmap[r.verdict] || "n/a";
    const row = document.createElement("button");
    row.className = "cand" + (r.passes ? "" : " blocked");
    row.innerHTML =
      `<span class="cand-rank">${i + 1}</span>
       <span class="cand-main">
         <span class="cand-name">${r.e3_name}<span class="cand-plus">+</span>${r.linker_name}</span>
         <span class="cand-bar"><span class="cand-bar-fill ${vc}" style="width:${r.score || 0}%"></span></span>
       </span>
       <span class="cand-meta">
         <span class="cand-score">${r.score == null ? "—" : r.score}</span>
         <span class="cand-mw">${r.mw} Da${r.passes ? "" : " · 🔒"}</span>
       </span>`;
    row.addEventListener("click", () => loadCombo(r.linker, r.e3));
    host.appendChild(row);
  });
}

function fillSelect(sel, items, label, allowNone) {
  sel.innerHTML = "";
  if (allowNone) {
    const o = document.createElement("option");
    o.value = ""; o.textContent = "None — raw developability";
    sel.appendChild(o);
  }
  items.forEach((it) => {
    const o = document.createElement("option");
    o.value = it.key; o.textContent = label(it);
    sel.appendChild(o);
  });
}

function fillLinkerSelect(items) {
  const sel = $("linker");
  sel.innerHTML = "";
  const groups = {};
  items.forEach((it) => { (groups[it.group] = groups[it.group] || []).push(it); });
  Object.entries(groups).forEach(([g, arr]) => {
    const og = document.createElement("optgroup");
    og.label = g;
    arr.forEach((it) => {
      const o = document.createElement("option");
      o.value = it.key; o.textContent = it.name;
      og.appendChild(o);
    });
    sel.appendChild(og);
  });
}

function metaFor(list, key) { return list.find((x) => x.key === key); }

function ordinal(n) {
  const s = ["th", "st", "nd", "rd"], v = n % 100;
  return n + (s[(v - 20) % 10] || s[v] || s[0]);
}

// ------------------------------------------------------------- update -------
let inflight = 0;
async function update() {
  const payload = {
    warhead: $("warhead").value,
    linker: $("linker").value,
    e3: $("e3").value,
    domain: $("domain").value || null,
  };

  // slot meta lines
  const w = metaFor(CAT.warheads, payload.warhead);
  const l = metaFor(CAT.linkers, payload.linker);
  const e = metaFor(CAT.e3, payload.e3);
  $("meta-warhead").textContent = w ? `${w.confidence} conf.` : "";
  $("meta-linker").textContent = l ? l.rigidity : "";
  $("meta-e3").textContent = e ? `${e.e3_class} · ${e.confidence} conf.` : "";

  $("loading").classList.add("show");
  const id = ++inflight;
  let data;
  try {
    data = await (await fetch("/api/assemble", {
      method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify(payload),
    })).json();
  } catch (err) {
    $("loading").classList.remove("show"); return;
  }
  if (id !== inflight) return; // a newer request superseded this one
  if (data.error) { $("loading").classList.remove("show"); return; }

  renderAll(data);
  $("loading").classList.remove("show");
}

async function uploadMol(ev) {
  const file = ev.target.files[0];
  if (!file) return;
  const text = await file.text();
  $("loading").classList.add("show");
  const id = ++inflight;
  let data;
  try {
    data = await (await fetch("/api/score_molblock", {
      method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ molblock: text, domain: $("domain").value || null }),
    })).json();
  } catch (e) { $("loading").classList.remove("show"); return; }
  ev.target.value = "";
  if (id !== inflight) return;
  if (data.error) { $("loading").classList.remove("show"); alert(data.error); return; }
  renderAll(data);
  $("loading").classList.remove("show");
}

function renderAll(data) {
  const uploaded = data.source === "uploaded";
  $("smiles-title").textContent = uploaded ? "Uploaded structure (SMILES)" : "Assembled SMILES";
  renderMolecule(data);
  renderMetrics(data);
  renderLinkerGeo(data.linker);
  renderTip(data.tip);
}

function renderLinkerGeo(lg) {
  const set = (id, txt, band) => { const el = $(id); el.textContent = txt; el.className = "geo-v " + (band || ""); };
  if (!lg || lg.error) { set("geo-len", "—"); set("geo-rot", "—"); set("geo-rig", "—"); return; }
  set("geo-len", lg.length_atoms, lg.length_band);
  set("geo-rot", lg.rot_bonds, lg.rot_band);
  set("geo-rig", lg.fsp3, lg.rigid_band);
}

function renderTip(tip) {
  const el = $("tip");
  if (tip) {
    $("tip-body").textContent = tip.replace(/^Tip:\s*/, "");
    el.classList.add("show");
  } else {
    el.classList.remove("show");
  }
}

// ------------------------------------------------------------- 3D -----------
function renderMolecule(data) {
  $("smiles").textContent = data.smiles;
  $("mw-pill").textContent = `${data.mw} Da · ${countAtoms(data.molblock)} atoms`;
  if (!data.molblock) { viewer.removeAllModels(); viewer.render(); return; }
  lastMolblock = data.molblock;
  viewer.removeAllModels();
  viewer.addModel(data.molblock, "sdf");
  applyStyle();
  viewer.zoomTo();
  viewer.render();
}

function applyStyle() {
  if (!lastMolblock) return;
  viewer.setStyle({}, styleSpec(curStyle));
  viewer.render();
}
function styleSpec(s) {
  if (s === "ball") return { stick: { radius: 0.14 }, sphere: { scale: 0.28 } };
  if (s === "line") return { line: { linewidth: 2 } };
  return { stick: { radius: 0.17 } };
}
function countAtoms(mb) {
  if (!mb) return 0;
  const line = mb.split("\n")[3] || "";
  const n = parseInt(line.slice(0, 3), 10);
  return isNaN(n) ? 0 : n;
}

// ------------------------------------------------------------- metrics ------
function renderMetrics(data) {
  // gauge
  const idx = data.index ?? 0;
  const fill = $("gauge-fill");
  fill.style.strokeDashoffset = GAUGE_LEN * (1 - idx / 100);
  fill.style.stroke = idx >= 70 ? COLORS.green : idx >= 45 ? COLORS.amber : COLORS.red;
  $("index-val").textContent = data.index ?? "—";

  // verdict
  const vmap = { Typical: "green", Borderline: "amber", Atypical: "red" };
  const vc = vmap[data.verdict] || "n/a";
  $("verdict-dot").style.background = COLORS[vc];
  $("verdict-txt").textContent = data.verdict;
  $("verdict-txt").style.color = COLORS[vc];

  // domain gate
  const gate = $("domain-gate");
  if (data.domain) {
    gate.classList.add("show");
    const d = data.domain;
    const chips = Object.entries(d.filter_results || {})
      .map(([k, ok]) => {
        const cls = ok === true ? "ok" : ok === false ? "no" : "";
        const mk = ok === true ? "✓" : ok === false ? "✕" : "?";
        return `<span class="gate-chip ${cls}">${mk} ${k}</span>`;
      }).join("");
    gate.innerHTML =
      `<div class="gate-head">Domain gate
         <span class="gate-badge ${d.passes ? "pass" : "block"}">${d.passes ? "PASS" : "BLOCKED"}</span>
         <span style="margin-left:auto;font-weight:600;color:var(--muted)">idx ${d.weighted_index}</span>
       </div>
       <div class="gate-filters">${chips || '<span class="gate-chip">no hard filters</span>'}</div>`;
  } else {
    gate.classList.remove("show");
  }

  // metric rows
  const host = $("metrics");
  host.innerHTML = "";
  data.metrics.forEach((m) => {
    const des = m.desirability;
    const pct = des == null ? 0 : Math.round(des * 100);
    const val = typeof m.value === "number"
      ? (Number.isInteger(m.value) ? m.value : m.value.toFixed(1))
      : (m.value ?? "n/a");
    const pctl = m.percentile;
    const pctlTxt = pctl == null ? "" :
      `<span class="m-pct" title="${ordinal(pctl)} percentile of PROTAC-DB (n=10,728)">${ordinal(pctl)} pctl</span>`;
    const row = document.createElement("div");
    row.className = "metric";
    row.innerHTML =
      `<span class="m-dot ${m.color}"></span>
       <div class="m-main">
         <div class="m-label">${m.label}</div>
         <div class="m-bar"><div class="m-bar-fill ${m.color}" style="width:${pct}%"></div></div>
       </div>
       <div class="m-val">
         <span class="m-num">${val}<span class="m-unit">${m.unit || ""}</span></span>
         ${pctlTxt}
       </div>`;
    host.appendChild(row);
  });
}

boot();
