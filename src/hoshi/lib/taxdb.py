# taxon_lineage.py
import taxopy
from collections import OrderedDict
from pathlib import Path
import os

# Global cache for TaxDb
_TAXDB_CACHE = {}


def get_taxdb(hoshi_path=None):
    """
    Return a cached TaxDb instance, optionally with a custom path.
    """
    global _TAXDB_CACHE
    HOSHI_DEFAULT = os.path.join(Path.home(), ".hoshi", "taxdb")
    # Use default path if not provided
    if hoshi_path is None:
        hoshi_path = HOSHI_DEFAULT

    # Normalize path
    hoshi_path = os.path.abspath(hoshi_path)

    if hoshi_path not in _TAXDB_CACHE:
        taxdb = taxopy.TaxDb(
            nodes_dmp=os.path.join(hoshi_path, "nodes.dmp"),
            names_dmp=os.path.join(hoshi_path, "names.dmp"),
        )
        _TAXDB_CACHE[hoshi_path] = taxdb

    return _TAXDB_CACHE[hoshi_path]


def get_filtered_lineage(taxid, taxdb=None, desired_ranks=None):
    """
    Return a filtered OrderedDict of ranks from species to kingdom for a taxon ID.
    """
    if taxdb is None:
        taxdb = get_taxdb()

    if desired_ranks is None:
        desired_ranks = [
            "species",
            "genus",
            "family",
            "order",
            "class",
            "phylum",
            "kingdom",
        ]

    taxon = taxopy.Taxon(taxid, taxdb)
    data = taxon.rank_name_dictionary
    return OrderedDict((k, data[k]) for k in desired_ranks if k in data)


def lineage_to_tsv(filtered_lineage, desired_ranks):
    """
    Format an OrderedDict lineage into a tab-delimited string based on desired_ranks.
    """
    return "\t".join(filtered_lineage.get(rank, "NA") for rank in desired_ranks)
