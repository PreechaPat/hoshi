"""Savont ingress — read Savont output directories into a ``SummarizedExperiment``.

Savont folder / species quirk
-----------------------------
Savont output has the same fixed directory structure (no sample-specific information),
so the sample name must be supplied explicitly via ``sample_names``.

Each Savont sample directory contains, after the initial (denoise/cluster)
run::

    feature-table.tsv    # ASV read counts            (used here)
    final_asvs.fasta     # representative ASV sequences
    final_clusters.tsv   # cluster membership

and, after the classify run, additional taxonomy outputs including::

    asv_mappings.tsv        # per-ASV taxonomy assignments  (used here)
    genus_abundance.tsv     # relative abundance at genus level
    species_abundance.tsv   # relative abundance at species level

Species resolution quirk
------------------------
1. The ASV mapping file can list multiple taxonomy hits per ASV. We resolve each
ASV to a **single species by taking the first hit** (``keep="first"``), then
aggregate read counts by ``tax_id`` to compute relative abundance.
"""

from __future__ import annotations

from pathlib import Path

from hoshi.lib.experiment import SummarizedExperiment
from hoshi.lib.savont_reader import SavontReader


def read_savont_abundance_into_summarizedexperiment(
    input_dir: str | Path,
    *,
    sample_name: str,
) -> SummarizedExperiment:
    """
    Load a single Savont output directory into a SummarizedExperiment.

    Input layout (Savont)
    ---------------------
    Savont output is 100% the same fixed directory structure for every sample so
    the sample name must be supplied explicitly via ``sample_name``.

    Each Savont sample directory contains, after the initial (denoise/cluster)
    run:
      - ``feature-table.tsv``   (ASV read counts — used here)
      - ``final_asvs.fasta``    (representative ASV sequences)
      - ``final_clusters.tsv``  (cluster membership)
    and, after the classify run, additional taxonomy outputs including:
      - ``asv_mappings.tsv``    (per-ASV taxonomy assignments — used here)
      - ``genus_abundance.tsv``  (relative abundance at the genus level)
      - ``species_abundance.tsv`` (relative abundance at the species level)

    Only ``feature-table.tsv`` and ``asv_mappings.tsv`` are consumed today. The
    other outputs (``final_asvs.fasta``, ``final_clusters.tsv``, and any further
    classify-stage files) are intentionally noted here because we will need them
    later for more detailed reports.

    This reads exactly one sample and returns a single-sample experiment whose
    feature index carries Savont's **raw** feature ids. Combining several
    samples into one experiment (and the sample-scoping that entails) is a
    separate, explicit step — see
    :func:`hoshi.lib.experiment.combine_experiments`.

    Parameters
    ----------
    input_dir : str or Path
        Path to a Savont output directory containing ``feature-table.tsv`` and
        ``asv_mappings.tsv``.

    sample_name : str
        Sample name for the resulting single-sample experiment.

    Returns
    -------
    SummarizedExperiment
        Single-sample container with:
        - assays["abundance"]: relative abundance matrix (features × 1)
        - assays["counts"]: estimated counts matrix (features × 1)
        - row_data: taxonomy annotations per feature (raw feature id index)
        - col_data: sample metadata (indexed by sample name)
        - metadata: {"source": "savont", "confidence": {...}}
    """
    return SavontReader(
        Path(input_dir), sample_name=sample_name
    ).to_summarized_experiment()
