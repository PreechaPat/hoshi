# Hoshi — Behavior Specification (SPECS)

> **Purpose.** This is the *behavior / requirements* spec for Hoshi, written to
> guide feature work — including AI-assisted development. Where [`SPEC.md`](SPEC.md)
> describes **how the system is built** (architecture, data model, layout), this
> document describes **what the system should do** (observable behavior,
> report contents, acceptance criteria).
>
> **Status legend**
> - ✅ **Implemented** — exists today, verified in code/tests.
> - 🚧 **Partial** — works but incomplete ("barebone").
> - 📋 **Planned** — agreed direction, not built.
> - ❓ **DECISION NEEDED** — open question; a human must decide before building.
>
> **How to use this with an AI.** Point the AI at a single requirement section,
> its acceptance criteria, and the relevant module. Ask it to (1) restate the
> criteria, (2) propose a plan, (3) implement, (4) prove each criterion with a
> test. Do not let it invent behavior for a ❓ item — surface the question first.

---

## 0. Scope & Vocabulary

- **Classifier** — an upstream 16S rRNA taxonomic tool (EMU, Savont, …).
- **SummarizedExperiment (SE)** — the internal standard: pure measurement data
  (abundance/counts assays + aligned taxonomy + sample metadata). No
  interpretation, no clinical fields. See `SPEC.md §1`.
- **Report** — the composite that pairs an SE with report-time data (confidence,
  clinical metadata) and renders to HTML/PDF. See `SPEC.md §3`.
- **QC metrics** — quality/run statistics *not* produced by the abundance step
  (read counts, quality scores, filtering stats, run identity). The primary
  driver of the "complete report" goal.
- **Complete report** — a report enriched beyond the current barebone templates
  with QC metrics and other-tool statistics, so a reader can judge the run, not
  just the taxonomy. This is the headline goal this SPECS tracks.

---

## 1. Product Goals

| # | Goal | Status |
|---|------|--------|
| G1 | Read heterogeneous classifier output into one SE. | ✅ (EMU, Savont) |
| G2 | Render single / multi / medical reports from an SE. | 🚧 barebone |
| G3 | **Enrich reports with QC + other-tool metrics ("complete report").** | 📋 this SPECS |
| G4 | Add classifiers without touching report code. | ✅ via `AbundanceReader` |
| G5 | Deterministic, testable report output. | 🚧 partial |

---

## 2. Current Behavior (baseline — verified in code)

These are the behaviors a "complete report" builds *on top of*. They are the
regression floor: changes must not break them.

### 2.1 `report-single` ✅ / 🚧
- **Input:** a classifier output directory (`--input-format savont` default, or
  `emu`).
- **Output:** responsive HTML at `-o` (default `report_single.html`); optional
  A4 PDF via `--pdf`/`--output-pdf` (WeasyPrint).
- **Contents today:** species table (top 5 / top 10), alpha diversity, an
  interactive taxonomy Sankey. Per-species confidence shown as a column + a
  headline percentage when the reader supplies it (Savont does, EMU does not).
- **Meta rows** (`unmapped`, `mapped_filtered`, `mapped_unclassified`) are
  filtered out of the organism display.

### 2.2 `report-multi` ✅ / 🚧
- **Input:** one or more classifier output directories.
- **Output:** one combined HTML, one tab per sample (species table + Sankey per
  tab). PDF supported but requires `-o/--output`.

### 2.3 `report-medical` ✅ / 🚧
- **Input:** a single-sample classifier output (EMU TSV *file* or Savont *dir*)
  plus an optional clinical metadata JSON (`-m`).
- **Output:** clinical-style HTML; optional PDF.
- **Clinical fields** come from the JSON; missing ones render as `N/A`, and an
  absent `authorized_by` renders as a blank signature line.
- **Pathogen classification** per organism from `--pathogen-sheet` CSV keyed by
  NCBI tax_id (commensal / potential / opportunistic / primary). Empty value
  disables the column.
- **Conclusion** taken from JSON, else auto-derived:
  no organisms → `not_detected`; any non-commensal → `pathogen_detected`;
  all commensal → `organism_detected`.
- **QC items:** if the JSON omits `qc_items`, defaults to two passing rows
  (`Read quality: pass`, `Read support: pass`). *(This is a placeholder — see §4.)*

### 2.4 `convert` ✅ / `enrich` ✅
- `convert`: EMU TSV → Kraken2 report (stdout or `-o`).
- `enrich`: add lineage columns to a `tax_id` table via the NCBI taxdb.

---

## 3. The "Complete Report" — Requirements (📋 primary work)

The core requirement: **surface QC and other-tool statistics in the report**, so
the report reflects run quality and not just taxonomy. This section defines *what*
the complete report should contain and the acceptance criteria; it deliberately
leaves *where the data comes from* and *how it is modeled* as ❓ items for you to
decide (you mentioned reading `SPEC.md` first).

### 3.1 QC Section (📋)

**Requirement.** Every report type gains a QC section that presents run/read
quality metrics sourced from tools upstream of (or alongside) the classifier.

**Candidate metrics** *(❓ confirm which apply to your pipeline)*:

| Metric | Example source | Notes |
|--------|----------------|-------|
| Total reads (raw) | sequencer / demux | per sample |
| Reads passing filter | QC tool (e.g. NanoPlot, fastp, Chopper) | count + % |
| Mean/median read length | QC tool | 16S full-length ≈ 1.5 kb |
| Mean quality (Q score) | QC tool | |
| Reads mapped / classified | classifier | already partly in SE meta rows |
| Reads unclassified / unmapped | classifier | already in SE meta rows |
| Negative-control signal | negative sample | contamination check |

> ❓ **DECISION NEEDED (D1): QC source of truth.** Where do QC metrics come from?
> Options: (a) a separate QC file per sample (JSON/TSV) passed via a new flag;
> (b) parsed from a known QC tool's native output (which tool?); (c) folded into
> the existing medical metadata JSON; (d) computed from data Hoshi already has
> (e.g. classifier meta rows). Likely a mix. Decide before modeling.

> ❓ **DECISION NEEDED (D2): where QC lives in the data model.** `SPEC.md §1.5`
> keeps the SE as pure measurement data, and `SPEC.md §4` proposes a planned
> `ReportContext` companion for non-measurement inputs. QC metrics are a natural
> fit for `ReportContext`. Confirm QC goes there (not into `SE.metadata`/`col_data`).

**Acceptance criteria (once D1/D2 are fixed):**
- [ ] Given a sample with QC data available, `report-single`/`-medical` renders a
      QC section listing the confirmed metrics with values and units.
- [ ] Given a sample with **no** QC data, the report renders without the QC
      section (or with an explicit "QC not available" note) — no crash, no fake
      values. *(Contrast with today's placeholder passing rows in §2.3.)*
- [ ] QC values are formatted deterministically (fixed rounding, thousands
      separators as agreed) so output is diff-stable for tests.
- [ ] A regression test asserts specific QC values appear in the generated HTML
      for a fixture sample.

### 3.2 Run / Read Statistics Panel (📋)
**Requirement.** A compact header/summary panel per sample: total reads, %
classified, number of taxa detected, dominant organism. Some of this is derivable
from the SE today (meta rows + abundance); the rest depends on D1.

**Acceptance criteria:**
- [ ] Panel shows: total reads (if available), % reads classified, distinct taxa
      count (excluding meta rows), top organism + its abundance.
- [ ] Values reconcile with the species table (top organism matches).

### 3.3 Other-tool integration (📋 — beyond QC)
> ❓ **DECISION NEEDED (D3): which other tools feed the complete report?** You
> mentioned "other tools, such as QC." List the concrete tools/outputs you want
> represented (e.g. host-depletion stats, primer-trimming stats, a specific QC
> reporter). Each becomes a small reader + a report section with its own criteria.

### 3.4 Negative control / contamination (📋, optional)
> ❓ **DECISION NEEDED (D4):** Should the complete report flag likely contaminants
> by comparing a sample against a negative control? `test_data` already contains
> `savont-out-negative01` / `negative01`, suggesting this is on your mind. If yes,
> define the comparison rule (e.g. subtract/flag taxa present in the negative
> above a threshold).

---

## 4. Known Placeholders to Replace

These exist today as stopgaps and should be resolved by the complete-report work:

- **Medical QC defaults (§2.3).** `report_to_medical_data` injects two passing QC
  rows when none are supplied. A complete report should derive real QC or clearly
  mark it absent — never emit fabricated "pass" rows. *(Ref: `command/report_medical.py`)*
- **Stale console scripts.** `pyproject.toml` declares `hoshi-report` and
  `hoshi-summary` pointing at non-existent modules. Remove or implement.

---

## 5. Cross-cutting Requirements

### 5.1 Determinism & Testability
- Report generation for a fixed input must produce stable, diff-able output
  (fixed dates injectable, sorted tables, fixed rounding).
- Each report section that carries data has at least one test asserting the data
  reaches the HTML (per `AGENTS.md`: "ensure generated HTML contains sample
  organism IDs").

### 5.2 Graceful degradation
- Missing optional inputs (QC file, metadata JSON, pathogen sheet) must never
  crash a report; they degrade to an omitted/`N/A` section with no fabricated
  values.

### 5.3 Separation of concerns (must hold)
- SE stays pure measurement data. Interpretation, QC, and clinical fields live in
  the report/context layer, not on the SE. *(This constrains how §3 is built —
  see D2.)*

### 5.4 PDF parity
- Anything added to HTML must render acceptably in the WeasyPrint PDF path, with
  the known exception that interactive Plotly Sankey degrades to a static/omitted
  panel in PDF.

---

## 6. Out of Scope (for the complete-report milestone)
- New classifier readers (Bracken/QIIME2/BIOM) — tracked in `SPEC.md`.
- Comparison / differential report between two experiments — `SPEC.md`.
- Formalizing the `IngressAdapter` protocol — `SPEC.md`.

---

## 7. Open Decisions Summary (fill these in)

| ID | Decision | Owner | Resolution |
|----|----------|-------|------------|
| D1 | QC metrics source of truth | Preecha | _TBD_ |
| D2 | QC data model location (confirm `ReportContext`) | Preecha | _TBD_ |
| D3 | Which other tools feed the complete report | Preecha | _TBD_ |
| D4 | Negative-control / contamination flagging | Preecha | _TBD_ |
| D5 | Which report adopts the complete-report sections first | Preecha | _TBD_ |

---

## 8. Suggested First Increment (proposal — needs your OK)

To make progress without waiting on every decision, a low-risk first slice:

1. Resolve **D1/D2** for the *single* most available QC source (likely the
   classifier meta rows already in the SE: reads mapped/classified/unmapped).
2. Add a **Run/Read Statistics panel (§3.2)** to `report-single` using only data
   Hoshi already has — no new inputs, no new file formats.
3. Add a regression test asserting the panel's numbers appear in the HTML and
   reconcile with the species table.
4. Only then design the external-QC-file path (§3.1) once D1/D3 are settled.

> This sequences the work so the first PR touches no ❓ items that require new
> external formats — it just surfaces data already in the SE.
