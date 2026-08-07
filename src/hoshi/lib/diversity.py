"""Diversity indices and summary statistics for microbiome abundance data."""

from __future__ import annotations

import math
from dataclasses import dataclass

import pandas as pd


@dataclass(frozen=True)
class DiversityStats:
    """Summary statistics for a single microbiome sample."""

    total_reads: int
    classified_reads: int
    unclassified_reads: int
    species_richness: int
    shannon_index: float
    simpson_index: float
    evenness: float  # Pielou's evenness: H / ln(S)


def compute_diversity(df: pd.DataFrame) -> DiversityStats:
    """
    Compute diversity statistics from an EMU abundance DataFrame.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with at least 'abundance' and 'estimated counts' columns,
        and a 'tax_id' column. Rows with tax_id in ('unmapped', 'mapped_filtered',
        'mapped_unclassified') are treated as non-species entries.

    Returns
    -------
    DiversityStats
        Computed diversity metrics.
    """
    # Separate classified species from metadata rows
    meta_ids = {"unmapped", "mapped_filtered", "mapped_unclassified"}
    is_species = ~df["tax_id"].astype(str).isin(meta_ids)

    species_df = df[is_species].copy()
    meta_df = df[~is_species].copy()

    # Read counts
    species_counts = pd.to_numeric(species_df["estimated counts"], errors="coerce").fillna(0)
    meta_counts = pd.to_numeric(meta_df["estimated counts"], errors="coerce").fillna(0)

    classified_reads = int(species_counts.sum())
    unclassified_reads = int(meta_counts.sum())
    total_reads = classified_reads + unclassified_reads

    # Species richness (number of species with non-zero abundance)
    nonzero = species_counts[species_counts > 0]
    species_richness = len(nonzero)

    # Relative proportions (for diversity indices, use only classified reads)
    if classified_reads > 0:
        proportions = nonzero / classified_reads
    else:
        proportions = pd.Series(dtype=float)

    # Shannon index: H = -sum(p_i * ln(p_i))
    shannon_index = -float((proportions * proportions.apply(math.log)).sum()) if len(proportions) > 0 else 0.0

    # Simpson index: 1 - sum(p_i^2) (inverse Simpson's diversity)
    simpson_index = 1.0 - float((proportions**2).sum()) if len(proportions) > 0 else 0.0

    # Pielou's evenness: J = H / ln(S)
    if species_richness > 1:
        evenness = shannon_index / math.log(species_richness)
    else:
        evenness = 0.0

    return DiversityStats(
        total_reads=total_reads,
        classified_reads=classified_reads,
        unclassified_reads=unclassified_reads,
        species_richness=species_richness,
        shannon_index=round(shannon_index, 4),
        simpson_index=round(simpson_index, 4),
        evenness=round(evenness, 4),
    )
