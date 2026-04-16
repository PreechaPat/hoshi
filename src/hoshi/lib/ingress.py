import os

import pandas as pd


def read_input_table(input_data, required_columns=None, sep=None):
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

    return df


def read_emu_abundance(input_data, *, sep: str | None = "\t", reorder: bool = False) -> pd.DataFrame:
    """
    Load an EMU rel-abundance table and normalise its column layout.
    """

    EMU_REQUIRED_COLUMNS = ("tax_id", "abundance")
    EMU_TAXONOMY_COLUMNS = (
        "superkingdom",
        "phylum",
        "class",
        "order",
        "family",
        "genus",
        "species",
    )
    EMU_OUTPUT_COLUMNS = (
        "tax_id",
        "abundance",
        "estimated counts",
        *EMU_TAXONOMY_COLUMNS,
    )
    EMU_NUMERIC_COLUMNS = ("abundance", "estimated counts")
    df = read_input_table(input_data, required_columns=list(EMU_REQUIRED_COLUMNS), sep=sep)

    # Only pick up those superkingdom -> species taxonomy
    for column in EMU_TAXONOMY_COLUMNS:
        if column not in df.columns:
            df[column] = pd.NA

    for column in EMU_NUMERIC_COLUMNS:
        if column in df.columns:
            df[column] = pd.to_numeric(df[column], errors="coerce")

    if reorder:
        if "estimated counts" in df.columns:
            df = df.rename(columns={"estimated counts": "estimated counts"})
            df["estimated counts"] = pd.to_numeric(df["estimated counts"], errors="coerce").astype(int)
        if "estimated counts" not in df.columns:
            df["estimated counts"] = pd.NA

        for column in EMU_OUTPUT_COLUMNS:
            if column not in df.columns:
                df[column] = pd.NA

        df = df.loc[:, list(EMU_OUTPUT_COLUMNS)]

    return df
