from __future__ import annotations

from pathlib import Path
from typing import NamedTuple

import pandas as pd

from hoshi.lib.ingress import read_emu_abundance

TAXONOMY_LEVELS = (
    "superkingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
)
VALUE_COLUMN_PRIORITY = ("estimated counts", "abundance")


class SankeyInput(NamedTuple):
    """Data required to render a Sankey plot."""

    frame: pd.DataFrame
    levels: list[str]
    value_column: str


def load_sankey_input(tsv_file: str | Path) -> SankeyInput:
    """Load and normalize an EMU abundance table ready for Sankey plotting."""

    dataframe = read_emu_abundance(str(tsv_file), reorder=True)
    return prepare_sankey_dataframe(dataframe)


def prepare_sankey_dataframe(df: pd.DataFrame) -> SankeyInput:
    """Normalize a dataframe of EMU abundances for Sankey plotting."""

    levels = _available_levels(df)
    value_column = _select_value_column(df)

    trimmed = df.loc[:, levels + [value_column]].copy()
    trimmed[value_column] = pd.to_numeric(trimmed[value_column], errors="coerce").fillna(0.0)
    trimmed = trimmed[trimmed[value_column] > 0]

    if trimmed.empty:
        raise ValueError("No taxonomic entries with positive abundance were found.")

    level_frame = trimmed[levels].fillna("").astype(str)
    trimmed.loc[:, levels] = level_frame
    trimmed = trimmed[~level_frame.eq("").all(axis=1)]

    if trimmed.empty:
        raise ValueError("All entries lack taxonomy assignments required for a Sankey plot.")

    return SankeyInput(frame=trimmed, levels=levels, value_column=value_column)


def _available_levels(df: pd.DataFrame) -> list[str]:
    levels = [level for level in TAXONOMY_LEVELS if level in df.columns]
    if len(levels) < 2:
        raise ValueError("Unable to derive two or more taxonomy levels from the EMU table.")
    return levels


def _select_value_column(df: pd.DataFrame) -> str:
    for column in VALUE_COLUMN_PRIORITY:
        if column in df.columns and df[column].notna().any():
            return column
    raise ValueError("Unable to locate a populated abundance column in the EMU table.")
