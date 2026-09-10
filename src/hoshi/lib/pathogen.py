"""Pathogen database — reads ``pathogen_sheet.csv`` into a queryable object.

The sheet is a curated CSV of known pathogens with NCBI ``taxid`` as the
primary key. Rather than passing raw dicts around, callers construct a
:class:`PathogenDB` once and query it by ``tax_id``:

    db = PathogenDB.from_csv("assets/pathogen_sheet.csv")
    db.classification("1496")        # -> "opportunistic"
    db.commensal_sites("817")        # -> "gut, stool"
    db.describe("1496")              # -> "opportunistic"
    db.describe("817")               # -> "commensal (gut, stool)"

The DB centralises how a sheet row is *encoded* into the single string the
report shows, so the presentation rule (e.g. appending commensal sites) lives
in one place instead of being scattered across callers/templates.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

# Columns the DB relies on. taxid + general_classification are required; the
# rest are optional and simply yield empty values when absent.
_REQUIRED_COLUMNS = ("taxid", "general_classification")


def _clean(value: object) -> str:
    """Normalise a cell to a stripped string ("" for NaN/None)."""
    if value is None or pd.isna(value):
        return ""
    return str(value).strip()


@dataclass(frozen=True)
class PathogenRecord:
    """A single pathogen sheet row, reduced to the fields we surface."""

    tax_id: str
    classification: str          # general_classification
    commensal_sites: str
    pathogenic_sites: str
    status: str                  # established / putative
    high_consequence: bool

    def describe(self) -> str:
        """Render the report-facing classification string.

        Commensal organisms carry their expected body sites in parentheses so a
        "commensal" call is not read as a bare verdict, e.g.
        ``"commensal (gut, stool)"``. Other classes render as the bare word.
        """
        if not self.classification:
            return ""
        if self.classification.lower() == "commensal" and self.commensal_sites:
            return f"{self.classification} ({self.commensal_sites})"
        return self.classification


class PathogenDB:
    """Queryable view over the pathogen sheet, keyed by ``tax_id``."""

    def __init__(self, records: dict[str, PathogenRecord]) -> None:
        self._records = records

    # ─── Construction ────────────────────────────────────────────────

    @classmethod
    def from_csv(cls, path: str | Path) -> PathogenDB:
        """Load a :class:`PathogenDB` from a pathogen sheet CSV."""
        df = load_pathogen_sheet(path)

        records: dict[str, PathogenRecord] = {}
        for _, row in df.iterrows():
            tax_id = _clean(row.get("taxid"))
            if not tax_id or tax_id in records:
                continue  # first record wins on duplicate tax_id
            records[tax_id] = PathogenRecord(
                tax_id=tax_id,
                classification=_clean(row.get("general_classification")),
                commensal_sites=_clean(row.get("commensal_sites")),
                pathogenic_sites=_clean(row.get("pathogenic_sites")),
                status=_clean(row.get("status")),
                high_consequence=_clean(row.get("high_consequence")).upper() == "TRUE",
            )
        return cls(records)

    # ─── Lookups ─────────────────────────────────────────────────────

    def __contains__(self, tax_id: object) -> bool:
        return str(tax_id) in self._records

    def __len__(self) -> int:
        return len(self._records)

    def record(self, tax_id: str) -> PathogenRecord | None:
        """Return the full record for ``tax_id`` (``None`` when unknown)."""
        return self._records.get(str(tax_id))

    def classification(self, tax_id: str) -> str | None:
        """Bare ``general_classification`` (``None`` when unknown)."""
        rec = self._records.get(str(tax_id))
        return rec.classification if rec else None

    def commensal_sites(self, tax_id: str) -> str | None:
        rec = self._records.get(str(tax_id))
        return rec.commensal_sites if rec else None

    def describe(self, tax_id: str) -> str | None:
        """Report-facing classification string (``None`` when unknown).

        See :meth:`PathogenRecord.describe` for the encoding rule.
        """
        rec = self._records.get(str(tax_id))
        return rec.describe() if rec else None

    def describe_lookup(self) -> dict[str, str]:
        """Flat ``{tax_id: describe()}`` map for all records.

        Convenient for passing to ``build_medical_report(pathogens=...)`` where
        the value is combined with the organism table at render time.
        """
        return {tax_id: rec.describe() for tax_id, rec in self._records.items()}


def load_pathogen_sheet(path: str | Path) -> pd.DataFrame:
    """Load and clean the pathogen sheet CSV.

    Returns the full DataFrame with ``taxid`` as a stripped string column
    (not the index) for easy joining against abundance tables.
    """
    path = Path(path)
    if not path.is_file():
        raise ValueError(f"Pathogen sheet not found: {path}")

    df = pd.read_csv(path, dtype=str)

    missing = [c for c in _REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise ValueError(
            f"Pathogen sheet missing required columns: {missing}. "
            f"Found: {list(df.columns)}"
        )

    # Normalise the taxid column — strip whitespace, drop empty/NaN rows.
    df["taxid"] = df["taxid"].str.strip()
    df = df[df["taxid"].notna() & (df["taxid"] != "")].copy()

    return df.reset_index(drop=True)
