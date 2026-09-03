"""Unified abundance-reader interface.

A *reader* wraps one classifier's output (a Savont directory, an EMU TSV, ...)
and knows how to turn it into a :class:`SummarizedExperiment`. Higher layers
(report commands) depend only on this interface, so switching input formats is
just choosing which concrete reader to instantiate:

    reader = build_reader(input_format, path, sample_name=name)
    se = reader.to_summarized_experiment()

Concrete readers
----------------
- :class:`hoshi.lib.savont_reader.SavontReader` — one Savont output directory.
- :class:`hoshi.lib.emu_reader.EmuReader`       — one or more EMU TSV files.

Both expose :meth:`AbundanceReader.to_summarized_experiment`, the single method
the report pipeline needs.
"""

from __future__ import annotations

from pathlib import Path
from typing import Protocol, runtime_checkable

from hoshi.lib.experiment import SummarizedExperiment

# Input formats understood by :func:`build_reader`. ``savont`` is the default
# across the CLI; ``emu`` remains supported for legacy inputs.
SUPPORTED_INPUT_FORMATS = ("savont", "emu")
DEFAULT_INPUT_FORMAT = "savont"


@runtime_checkable
class AbundanceReader(Protocol):
    """Common interface for classifier-output readers.

    Any reader can be handed to the report pipeline as long as it can produce a
    :class:`SummarizedExperiment`.
    """

    def to_summarized_experiment(self) -> SummarizedExperiment:
        """Read the wrapped output into a :class:`SummarizedExperiment`."""
        ...


def build_reader(
    input_format: str,
    input_path: str | Path,
    *,
    sample_name: str | None = None,
) -> AbundanceReader:
    """Instantiate the reader for ``input_format``.

    Parameters
    ----------
    input_format : str
        One of :data:`SUPPORTED_INPUT_FORMATS` (``"savont"`` or ``"emu"``).
    input_path : str or Path
        For ``savont``: a Savont output directory. For ``emu``: an EMU
        ``*_rel-abundance.tsv`` file.
    sample_name : str, optional
        Sample name for the resulting experiment. Defaults to the directory
        name (Savont) or the filename-derived stem (EMU).

    Returns
    -------
    AbundanceReader
        A concrete reader ready for :meth:`to_summarized_experiment`.

    Raises
    ------
    ValueError
        If ``input_format`` is not supported.
    """
    # Imported lazily to keep heavy pandas-backed modules off callers that only
    # need the interface / constants.
    if input_format == "savont":
        from hoshi.lib.savont_reader import SavontReader  # noqa: PLC0415

        return SavontReader(input_path, sample_name=sample_name)

    if input_format == "emu":
        from hoshi.lib.emu_reader import EmuReader  # noqa: PLC0415

        return EmuReader(input_path, sample_name=sample_name)

    raise ValueError(
        f"Unsupported input format '{input_format}'. "
        f"Supported: {', '.join(SUPPORTED_INPUT_FORMATS)}"
    )
