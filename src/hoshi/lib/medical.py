"""Medical-report intermediate built on the :class:`Report` composite.

Flow (Option A — everything clinical lives in ``Report.metadata``):

    inputs (classifier output + clinical params)
        -> build_medical_report(...) -> Report        # intermediate
        -> report_to_medical_data(report) -> dict      # flattened for Jinja2
        -> template.render(**data)                     # HTML

Keeping the clinical envelope in ``metadata`` (rather than typed fields) means
the whole ``Report`` can be dumped to JSON later without a bespoke schema, and
per-organism signals such as ``pathogenic`` can be layered in as a plain
``{tax_id: bool}`` dict that is combined with the experiment-derived organism
table at render time.
"""

from __future__ import annotations

from typing import Any

import pandas as pd

from hoshi.lib.experiment import SummarizedExperiment
from hoshi.lib.report import Report

# tax_id values EMU uses for control/meta rows that must never be organisms.
_META_TAX_IDS = {"unmapped", "mapped_filtered", "mapped_unclassified"}

# Clinical fields the medical template consumes, with render-time defaults.
# Required fields have a value of ``_REQUIRED`` sentinel and are validated.
_REQUIRED = object()

_CLINICAL_FIELDS: dict[str, Any] = {
    "report_id": _REQUIRED,
    "patient_id": _REQUIRED,
    "specimen_id": _REQUIRED,
    "specimen_type": _REQUIRED,
    "collection_date": _REQUIRED,
    "conclusion": "not_detected",
    "status": "Final",
    "lab_name": None,
    "reference_db": None,
    "method": None,
    "authorized_by": None,
}


def build_medical_report(
    experiment: SummarizedExperiment,
    *,
    clinical: dict[str, Any] | None = None,
    pathogens: dict[str, bool] | None = None,
    qc_items: list[dict[str, str]] | None = None,
    confidence: dict[str, float] | None = None,
) -> Report:
    """Build the ``Report`` intermediate for a medical report.

    Parameters
    ----------
    experiment : SummarizedExperiment
        The abundance/taxonomy container for the sample. Organisms are derived
        from this at render time (kept pure — never mutated here).
    clinical : dict, optional
        Clinical envelope (report_id, patient_id, specimen_*, conclusion, ...).
        Stored verbatim in ``Report.metadata``.
    pathogens : dict[str, bool], optional
        Per-``tax_id`` pathogenicity flags, combined with the organism table at
        render time. Stored under ``metadata["pathogens"]``.
    qc_items : list[dict], optional
        QC rows (each ``{"name": ..., "status": "pass"|"fail"}``). Stored under
        ``metadata["qc_items"]``.
    confidence : dict[str, float], optional
        Per-``tax_id`` species-calling confidence (percent). Populates the
        typed ``Report.confidence`` field and surfaces as organism ``identity``.
    """
    metadata: dict[str, Any] = dict(clinical or {})
    if pathogens is not None:
        metadata["pathogens"] = dict(pathogens)
    if qc_items is not None:
        metadata["qc_items"] = list(qc_items)

    report = Report(experiment=experiment, metadata=metadata)
    if confidence:
        report = report.with_confidence(confidence)
    return report


def _organisms_from_report(report: Report, *, top: int | None = None) -> list[dict]:
    """Derive the organism table from the experiment, combining pathogens.

    Organisms are sorted by abundance descending. ``identity`` comes from
    ``Report.confidence`` and ``pathogenic`` from ``metadata["pathogens"]`` —
    both combined here at render time, keyed by ``tax_id``.
    """
    df = report.experiment.to_dataframe()

    if "tax_id" in df.columns:
        df = df[~df["tax_id"].astype(str).isin(_META_TAX_IDS)].copy()

    df = df[pd.to_numeric(df["abundance"], errors="coerce").notna()].copy()
    df = df.sort_values("abundance", ascending=False)
    if top is not None:
        df = df.head(top)

    pathogens = report.metadata.get("pathogens", {})
    confidence = report.confidence

    organisms: list[dict] = []
    for _, row in df.iterrows():
        tax_id = str(row["tax_id"])

        name = row.get("species", "")
        if not name or pd.isna(name):
            name = row.get("genus", "Unknown")
            if pd.isna(name):
                name = "Unknown"

        identity = confidence.get(tax_id)
        organisms.append(
            {
                "name": str(name).strip(),
                "abundance": round(float(row["abundance"]) * 100, 2),
                "identity": round(identity, 1) if identity is not None else None,
                "pathogenic": pathogens.get(tax_id),
            }
        )

    return organisms


def report_to_medical_data(report: Report, *, top: int | None = None) -> dict:
    """Flatten a medical ``Report`` into the dict the template consumes.

    Required clinical fields (report_id, patient_id, specimen_id,
    specimen_type, collection_date) must be present in ``report.metadata`` —
    this mirrors the validation the JSON path performs.
    """
    meta = report.metadata

    missing = [
        key
        for key, default in _CLINICAL_FIELDS.items()
        if default is _REQUIRED and key not in meta
    ]
    if missing:
        raise ValueError(
            f"Missing required clinical fields in Report.metadata: {', '.join(missing)}"
        )

    data: dict[str, Any] = {}
    for key, default in _CLINICAL_FIELDS.items():
        data[key] = meta[key] if default is _REQUIRED else meta.get(key, default)

    data["organisms"] = _organisms_from_report(report, top=top)

    qc_items = meta.get("qc_items")
    if qc_items is not None:
        data["qc_items"] = list(qc_items)

    return data
