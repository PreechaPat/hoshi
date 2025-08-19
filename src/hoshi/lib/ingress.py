import pandas as pd
import os


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
            df = pd.read_csv(
                input_data, sep=sep, engine="python"
            )  # auto-detect by default
        except Exception as e:
            raise ValueError(f"Failed to read file '{input_data}': {e}")
    else:
        raise ValueError("input_data must be a file path or a pandas DataFrame.")

    if required_columns:
        missing = [col for col in required_columns if col not in df.columns]
        if missing:
            raise ValueError(f"Missing required columns: {missing}")

    return df
