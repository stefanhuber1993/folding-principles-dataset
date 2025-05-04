import numpy as np
import pandas as pd
from dataset.geometry import get_ca, get_cb, identify_strands

def detect_hairpins_and_chirality(model, dssp_df):
    """
    Detect β-hairpins in a DSSP-annotated structure and compute their handedness.

    Parameters:
        model : Bio.PDB.Model.Model
            The model object of the PDB structure (e.g., structure[0])
        dssp_df : pd.DataFrame
            DataFrame output from DSSP containing simplified secondary structure

    Returns:
        pd.DataFrame with β-hairpin metadata including chirality
    """
    pdb_id = getattr(getattr(model, "parent", None), "id", "UNKNOWN")
    hairpins = []

    for chain_id in dssp_df["Chain"].unique():
        chain_df = dssp_df[dssp_df["Chain"] == chain_id].reset_index(drop=True)
        strands = identify_strands(chain_df)

        for i in range(len(strands) - 1):
            s1_start, s1_end = strands[i]
            s2_start, s2_end = strands[i + 1]

            # Check that only coil residues exist in the loop region
            if np.any(chain_df.loc[s1_end + 1 : s2_start - 1, "SimpleSS"] != "C"):
                continue

            ca1 = get_ca(model, chain_id, chain_df.loc[s1_start, "ResNum"])
            ca2 = get_ca(model, chain_id, chain_df.loc[s1_end,   "ResNum"])
            ca3 = get_ca(model, chain_id, chain_df.loc[s2_start, "ResNum"])
            ca4 = get_ca(model, chain_id, chain_df.loc[s2_end,   "ResNum"])
            cb_pre = get_cb(model, chain_id, chain_df.loc[s1_end, "ResNum"])

            if None in (ca1, ca2, ca3, ca4, cb_pre):
                continue

            u = ca2 - ca1
            mid1 = (ca1 + ca2) / 2.0
            mid2 = (ca3 + ca4) / 2.0
            v = mid2 - mid1
            n = cb_pre - ca2

            sign = np.dot(np.cross(u.get_array(), v.get_array()), n.get_array())
            handedness = "R" if sign > 0 else "L"

            hairpins.append({
                "PDB": pdb_id,
                "Chain": chain_id,
                "strand1_start": chain_df.loc[s1_start, "ResNum"],
                "strand1_end":   chain_df.loc[s1_end,   "ResNum"],
                "strand2_start": chain_df.loc[s2_start, "ResNum"],
                "strand2_end":   chain_df.loc[s2_end,   "ResNum"],
                "loop_len":      s2_start - s1_end - 1,
                "handedness":    handedness
            })

    return pd.DataFrame(hairpins)
