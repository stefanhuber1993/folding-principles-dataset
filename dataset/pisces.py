import pandas as pd
import requests
from io import StringIO

DEFAULT_PISCES_URL = (
    "http://dunbrack.fccc.edu/pisces/download/"
    "cullpdb_pc25.0_res0.0-2.5_len40-10000_R0.3_Xray_d2025_02_19_chains11652"
)

def fetch_pisces_table(url: str = DEFAULT_PISCES_URL) -> pd.DataFrame:
    """Fetch and parse PISCES structure list into a pandas DataFrame."""
    r = requests.get(url, timeout=10)
    r.raise_for_status()
    lines = r.text.splitlines()

    # Extract column headers from the first data line
    for line in lines:
        if not line.startswith("#"):
            headers = line.split()
            break

    df = pd.read_csv(StringIO("\n".join(lines)), sep=r"\s+", comment="#", names=headers, header=None)
    df = df.iloc[1:].reset_index(drop=True)

    # Show available columns for debugging
    print("Detected PISCES columns:", df.columns.tolist())

    # Normalize expected columns
    colmap = {
        "pdb": "pdb",
        "PDB": "pdb",
        headers[0]: "pdb",  # force first column to be PDB code
    }
    df = df.rename(columns=colmap)

    df["resol"] = pd.to_numeric(df["resol"], errors="coerce")
    df["freerfac"] = pd.to_numeric(df["freerfac"], errors="coerce")
    df["len"] = pd.to_numeric(df["len"], errors="coerce")

    return df



def get_filtered_pdb_codes(
    df: pd.DataFrame,
    max_resolution: float = 2.5,
    max_rfree: float = 0.30,
    min_length: int = 40,
    max_length: int = 10000
) -> list[str]:
    """
    Filter the PISCES DataFrame and return a list of unique 4-letter PDB codes.

    Parameters:
        df : DataFrame returned by fetch_pisces_table()
        max_resolution : float (Å)
        max_rfree : float
        min_length : int
        max_length : int

    Returns:
        List of uppercase 4-letter PDB codes (chain IDs removed).
    """
    filtered = df[
        (df["resol"] <= max_resolution) &
        (df["freerfac"] <= max_rfree) &
        (df["len"] >= min_length) &
        (df["len"] <= max_length)
    ].copy()

    # Extract 4-letter PDB codes only
    filtered["pdb"] = filtered["pdb"].str[:4].str.upper()

    return filtered["pdb"].drop_duplicates().tolist()