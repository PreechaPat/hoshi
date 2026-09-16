"""EmuReader — an interface to a single EMU output directory.

The Savont/EMU structural difference (see :mod:`hoshi.lib.ingress_emu`): EMU
emits **one TSV file per sample** and encodes the sample name in the filename,
whereas Savont uses a fixed per-sample *directory*. To keep a consistent reader
interface, ``EmuReader`` also wraps a **directory** — an EMU output folder that
contains one ``*_rel-abundance.tsv`` file (the abundance table consumed today)
and may *optionally* contain other EMU outputs (none are emitted by default),
e.g.::

    B1_rel-abundance.tsv                      (used, required)
    B1_read-assignment-distributions.tsv      (optional, not yet used)
    B1_emu_alignments.bam                      (optional, not yet used)
    B1_unclassified_mapped.fastq.gz            (optional, not yet used)
    B1_unmapped.fastq.gz                        (optional, not yet used)

The sample name defaults to the EMU filename prefix (see
:func:`_derive_sample_name`), which is where EMU carries sample identity.

This is the EMU counterpart to :class:`hoshi.lib.savont_reader.SavontReader`.
Both take a directory and expose :meth:`to_summarized_experiment`, so the report
pipeline can consume either through the
:class:`hoshi.lib.reader.AbundanceReader` interface.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property
from pathlib import Path

from hoshi.lib.experiment import SummarizedExperiment

# EMU output filename suffixes, keyed by role. Only ``abundance`` is required
# and consumed today; the rest are optional (EMU emits none of them by default)
# and are discovered so the reader owns the full folder layout and can surface
# them later without changing its interface.
_EMU_SUFFIXES = {
    "abundance": "_rel-abundance.tsv",
    "read_assignments": "_read-assignment-distributions.tsv",
    "unclassified": "_unclassified_mapped.fastq.gz",
    "unmapped": "_unmapped.fastq.gz",
}


def _derive_sample_name(path: Path) -> str:
    """Derive a sample name from an EMU output filename.

    EMU encodes the sample name as the filename prefix (e.g.
    ``sample01_rel-abundance.tsv`` -> ``sample01``,
    ``barcode11.fastq_rel-abundance.tsv`` -> ``barcode11``). This naming rule is
    EMU-specific, so it lives with :class:`EmuReader` rather than in the shared
    reader interface.
    """
    stem = path.stem
    # Strip common EMU suffixes.
    for suffix in ("_rel-abundance", ".fastq_rel-abundance", "_rel_abundance"):
        if stem.endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    # Also handle patterns like "barcode11.fastq_rel-abundance".
    if stem.endswith(".fastq"):
        stem = stem[: -len(".fastq")]
    return stem


@dataclass(frozen=True)
class EmuReader:
    """Interface to one EMU output directory.

    Parameters
    ----------
    directory : str or Path
        Path to an EMU output directory containing a ``*_rel-abundance.tsv``
        file (and optionally other EMU outputs).
    sample_name : str, optional
        Sample name used when producing a :class:`SummarizedExperiment`.
        Defaults to the EMU filename-derived stem (see ``_derive_sample_name``).

    Raises
    ------
    ValueError
        If ``directory`` is not a directory or contains no
        ``*_rel-abundance.tsv`` file.
    """

    directory: Path
    sample_name: str

    def __init__(self, directory: str | Path, *, sample_name: str | None = None) -> None:
        path = Path(directory)
        if not path.is_dir():
            raise ValueError(f"EMU input must be a directory: {path}")
        abundance = _find_emu_file(path, _EMU_SUFFIXES["abundance"])
        if abundance is None:
            raise ValueError(
                f"No *{_EMU_SUFFIXES['abundance']} file found in {path}. "
                "Is this an EMU output directory?"
            )
        object.__setattr__(self, "directory", path)
        object.__setattr__(
            self, "sample_name", sample_name or _derive_sample_name(abundance)
        )

    # ─── Fixed-layout file discovery ─────────────────────────────────

    @cached_property
    def files(self) -> dict[str, Path | None]:
        """Discovered EMU output files by role.

        ``abundance`` is guaranteed present (enforced in ``__init__``); every
        other role is optional and is ``None`` when absent (EMU emits none of
        them by default).
        """
        return {
            role: _find_emu_file(self.directory, suffix)
            for role, suffix in _EMU_SUFFIXES.items()
        }

    @property
    def abundance_path(self) -> Path:
        """Path to the ``*_rel-abundance.tsv`` file (guaranteed present)."""
        path = self.files["abundance"]
        assert path is not None  # enforced in __init__
        return path

    # ─── Unified reader interface ────────────────────────────────────

    def to_summarized_experiment(self) -> SummarizedExperiment:
        """Build a single-sample :class:`SummarizedExperiment` for this folder.

        EMU has no per-call identity signal, so (unlike Savont) the result
        carries no ``species_confidence`` in its metadata.

        Returns
        -------
        SummarizedExperiment
            Container with:
            - assays["abundance"]: relative abundance matrix (features × 1)
            - assays["counts"]: estimated counts matrix (features × 1)
            - row_data: taxonomy annotations per feature (indexed by tax_id)
            - col_data: sample metadata (indexed by sample name)
            - metadata: {"source": "emu"}
        """
        # Imported lazily to keep heavy pandas-backed ingress off the interface.
        # Import via the ``ingress`` hub (not ``ingress_emu`` directly) so the
        # ingress<->ingress_emu re-export cycle is resolved in the right order.
        from hoshi.lib.ingress import (  # noqa: PLC0415
            read_emu_abundance_into_summarizedexperiment,
        )

        return read_emu_abundance_into_summarizedexperiment(
            self.abundance_path, sample_name=self.sample_name
        )


def _find_emu_file(directory: Path, suffix: str) -> Path | None:
    """Return the first file in ``directory`` whose name ends with ``suffix``."""
    for f in sorted(directory.iterdir()):
        if f.is_file() and f.name.endswith(suffix):
            return f
    return None
