"""Shared ingress utilities plus re-exports of source-specific readers.

The per-source ingress logic lives in dedicated modules so each can document its
own folder / species quirks without one file growing unwieldy:

- :mod:`hoshi.lib.ingress.emu`    — EMU rel-abundance TSVs (sample name is the
  filename prefix; no fixed directory layout).
- :mod:`hoshi.lib.ingress.savont` — Savont output directories (100% identical
  fixed directory structure; sample name supplied explicitly).

The EMU/Savont readers are re-exported from here so existing
``from hoshi.lib.ingress import ...`` imports keep working. This module itself
only owns :func:`read_input_table`, the generic delimited-table loader shared by
the source-specific readers.
"""

from __future__ import annotations

import os

import pandas as pd


def read_input_table(input_data, required_columns=None, sep=None, keep_only_required=False):
    """
    Read a TSV/CSV file or a pandas DataFrame and validate required columns.

    Parameters
    ----------
    input_data : str or pd.DataFrame
        Input file path or already-loaded DataFrame.

    required_columns : list of str, optional
        List of column names that must exist in the table.
        If any are missing, raises ValueError.

    sep : str, optional
        Delimiter for file reading (e.g., '\t' for TSV). If None, will auto-detect.

    keep_only_required : bool, default False
        When True (and ``required_columns`` is given), drop every column that is
        not in ``required_columns``, returning the required columns in the order
        they are listed. Lets callers hand off column pruning to the loader
        instead of slicing the frame themselves.

    Returns
    -------
    pd.DataFrame
        The loaded and validated DataFrame.

    Raises
    ------
    ValueError
        If required columns are missing or file cannot be parsed.
    """
    if isinstance(input_data, pd.DataFrame):
        df = input_data.copy()
    elif isinstance(input_data, str):
        if not os.path.isfile(input_data):
            raise ValueError(f"File not found: {input_data}")
        try:
            df = pd.read_csv(input_data, sep=sep, engine="python")  # auto-detect by default
        except Exception as e:
            raise ValueError(f"Failed to read file '{input_data}': {e}")
    else:
        raise ValueError("input_data must be a file path or a pandas DataFrame.")

    if required_columns:
        missing = [col for col in required_columns if col not in df.columns]
        if missing:
            raise ValueError(f"Missing required columns: {missing}")
        if keep_only_required:
            df = df.loc[:, list(required_columns)]

    return df


# ── Backward-compatible re-exports ───────────────────────────────────────────
# Imported at the bottom to avoid a circular import: ``ingress.emu`` imports
# ``read_input_table`` from this module, which must be defined first.
from hoshi.lib.ingress.emu import (  # noqa: E402
    read_emu_abundance,
    read_emu_abundance_into_summarizedexperiment,
    _derive_sample_name,
)
from hoshi.lib.ingress.savont import (  # noqa: E402
    read_savont_abundance,
    read_savont_abundance_into_summarizedexperiment,
)

__all__ = [
    "read_input_table",
    "read_emu_abundance",
    "read_emu_abundance_into_summarizedexperiment",
    "read_savont_abundance",
    "read_savont_abundance_into_summarizedexperiment",
    "_derive_sample_name",
]
