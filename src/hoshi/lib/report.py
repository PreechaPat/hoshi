"""Composite report object for 16S analysis.

``Report`` is a *composite* over :class:`SummarizedExperiment`: it holds a
reference to the underlying experiment (the pure abundance/taxonomy container)
and layers report-time analysis data on top of it — currently species-calling
confidence, plus a free-form metadata dict for anything else we need later.

The name ``Report`` is intentionally generic for now and may be renamed later.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
from typing import Any

from hoshi.lib.experiment import SummarizedExperiment


@dataclass()
class Report:
    """A 16S analysis report composed over a :class:`SummarizedExperiment`.

    Attributes
    ----------
    experiment : SummarizedExperiment
        The underlying abundance/taxonomy container. Kept pure — this object
        never mutates it.

    confidence : dict[str, float]
        Per-OTU calling confidence, keyed by the OTU/ASV feature id, expressed
        as a percent (0–100). For Savont this is each ASV's
        ``alignment_identity``. Sources without a per-call identity signal (e.g.
        EMU) leave this empty. Reports are rendered at species/genus level, so
        consumers aggregate these per-OTU values to the display taxon with
        ``max()`` (see ``aggregate_to_species``).

    metadata : dict[str, Any]
        Free-form report-level metadata (kept "just in case" for values that do
        not yet warrant a dedicated field).
    """

    experiment: SummarizedExperiment
    confidence: dict[str, float] = field(default_factory=dict)
    metadata: dict[str, Any] = field(default_factory=dict)

    # ─── Construction ────────────────────────────────────────────────

    @classmethod
    def from_experiment(
        cls,
        experiment: SummarizedExperiment,
        *,
        metadata: dict[str, Any] | None = None,
    ) -> Report:
        """Build a ``Report`` from an experiment, auto-extracting confidence.

        Per-OTU calling confidence is read from
        ``experiment.metadata["confidence"]`` when present (Savont populates it;
        EMU does not). For a single-sample experiment the sample's
        ``{feature_id: pct}`` map is used; multi-sample experiments start with no
        headline confidence (each sample is reported per-column downstream).

        This is the shared entry point for all report commands so they layer
        report-time data on the experiment the same way.
        """
        per_sample = experiment.metadata.get("confidence")
        confidence: dict[str, float] = {}
        if per_sample and experiment.n_samples == 1:
            sample = str(experiment.sample_ids[0])
            confidence = dict(per_sample.get(sample, {}))

        return cls(
            experiment=experiment,
            confidence=confidence,
            metadata=dict(metadata or {}),
        )

    # ─── Convenience ─────────────────────────────────────────────────

    @property
    def confidence_pct(self) -> float | None:
        """Experiment-level calling confidence, as a percent.

        Reports render at species/genus level, so this collapses the per-OTU
        confidence values to a single headline number using ``max()``.

        TODO: this is a placeholder aggregation. Revisit to use an
        abundance-weighted, species-level confidence rather than a plain max.
        """
        if not self.confidence:
            return None
        return max(self.confidence.values())

    # ─── Immutable update helpers ────────────────────────────────────

    def with_confidence(self, confidence: dict[str, float]) -> Report:
        """Return a copy with ``confidence`` replaced (no mutation)."""
        return replace(self, confidence=dict(confidence))

    def with_metadata(self, **updates: Any) -> Report:
        """Return a copy with ``metadata`` merged with ``updates`` (no mutation)."""
        merged = {**self.metadata, **updates}
        return replace(self, metadata=merged)

    # ─── Representation ──────────────────────────────────────────────

    def __repr__(self) -> str:
        pct = self.confidence_pct
        pct_str = f"{pct:.2f}%" if pct is not None else "n/a"
        return (
            f"Report("
            f"experiment={self.experiment!r}, "
            f"n_confidence={len(self.confidence)}, "
            f"confidence_pct={pct_str})"
        )
