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
        Species-calling confidence, keyed by ``tax_id``, expressed as a percent
        (0–100). For Savont this is the max ``alignment_identity`` across the
        ASVs that resolve to each species. Sources without a per-call identity
        signal (e.g. EMU) leave this empty.

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

        Species-calling confidence is read from
        ``experiment.metadata["species_confidence"]`` when present (Savont
        populates it; EMU does not). For a single-sample experiment the sample's
        ``{tax_id: pct}`` map is used; multi-sample experiments start with no
        headline confidence (each sample is reported per-column downstream).

        This is the shared entry point for all report commands so they layer
        report-time data on the experiment the same way.
        """
        per_sample = experiment.metadata.get("species_confidence")
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
        """Experiment-level species-calling confidence, as a percent.

        Computed as the abundance-weighted mean of the per-species confidence
        values — species that make up more of the sample count more toward the
        headline number. Falls back to a plain mean when abundance weights are
        unavailable, and returns ``None`` when there is no confidence data.
        """
        if not self.confidence:
            return None

        weights = self._abundance_weights()
        if weights:
            total_weight = sum(weights.get(tax_id, 0.0) for tax_id in self.confidence)
            if total_weight > 0:
                weighted = sum(
                    value * weights.get(tax_id, 0.0)
                    for tax_id, value in self.confidence.items()
                )
                return weighted / total_weight

        # Fallback: unweighted mean.
        return sum(self.confidence.values()) / len(self.confidence)

    def _abundance_weights(self) -> dict[str, float]:
        """Per-tax_id abundance weight summed across samples (empty if absent)."""
        assays = self.experiment.assays
        if "abundance" not in assays:
            return {}
        abundance = assays["abundance"]
        # Sum across samples so multi-sample reports weight by total abundance.
        per_feature = abundance.sum(axis=1)
        return {str(tax_id): float(value) for tax_id, value in per_feature.items()}

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
