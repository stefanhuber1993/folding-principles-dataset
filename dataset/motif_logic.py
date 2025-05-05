import numpy as np
from Bio.PDB import Vector
from .geometry import get_ca_from_residue, get_cb_from_residue
import pandas as pd

def identify_strands(dssp_df, min_len=2):
    """
    Identify continuous β-strand segments (E) in a DSSP-annotated DataFrame.

    Parameters:
        dssp_df : pd.DataFrame
        min_len : int — minimum length of strand segments

    Returns:
        List of (start_idx, end_idx) tuples for each β-strand.
    """
    strands = []
    in_seg = False
    for i, row in dssp_df.iterrows():
        if row["SimpleSS"] == "E":
            if not in_seg:
                seg_start = i
                in_seg = True
        else:
            if in_seg and (i - seg_start) >= min_len:
                strands.append((seg_start, i - 1))
                in_seg = False
    if in_seg and (len(dssp_df) - seg_start) >= min_len:
        strands.append((seg_start, len(dssp_df) - 1))
    return strands


def compute_handedness_from_strand_vectors(first_strand_axis, mid1, mid2, sidechain_vector):
    """
    Compute β-hairpin handedness using strand vectors and a sidechain vector.

    Parameters:
        first_strand_axis : Vector from start to end of first strand
        mid1, mid2        : Vectors for strand centers
        sidechain_vector  : CA→CB of residue preceding the loop

    Returns:
        handedness ("L"/"R") and normalized scalar triple product
    """
    u = first_strand_axis.get_array()
    v = (mid2 - mid1).get_array()
    n = sidechain_vector.get_array()
    numerator = np.dot(np.cross(u, v), n)
    denominator = np.linalg.norm(u) * np.linalg.norm(v) * np.linalg.norm(n)
    magnitude = numerator / denominator if denominator > 1e-6 else 0.0
    handedness = "R" if magnitude > 0 else "L"
    return handedness, magnitude


def detect_hairpins_and_chirality(model, dssp_df):
    """
    Detect β-hairpins and compute handedness using strand-axis and center geometry.

    This is the "standard" method using strand axis and strand midpoints.

    Parameters:
        model    : Bio.PDB Model object
        dssp_df  : DSSP output DataFrame with simplified secondary structure

    Returns:
        pd.DataFrame with per-hairpin metadata and chirality assignments
    """
    pdb_id = getattr(getattr(model, "parent", None), "id", "UNKNOWN")
    hairpins = []

    for chain_id in dssp_df["Chain"].unique():
        chain_df = dssp_df[dssp_df["Chain"] == chain_id].reset_index(drop=True)
        strands = identify_strands(chain_df)
        chain_sequence = ''.join(chain_df["AA"].tolist())

        for i in range(len(strands) - 1):
            s1_start, s1_end = strands[i]
            s2_start, s2_end = strands[i + 1]

            # Require uninterrupted coil between strands
            if np.any(chain_df.loc[s1_end + 1 : s2_start - 1, "SimpleSS"] != "C"):
                continue

            try:
                res_s1_start = model[chain_id][(" ", chain_df.loc[s1_start, "ResNum"], " ")]
                res_s1_end   = model[chain_id][(" ", chain_df.loc[s1_end,   "ResNum"], " ")]
                res_s2_start = model[chain_id][(" ", chain_df.loc[s2_start, "ResNum"], " ")]
                res_s2_end   = model[chain_id][(" ", chain_df.loc[s2_end,   "ResNum"], " ")]
            except KeyError:
                continue

            # Fetch geometry
            ca1 = get_ca_from_residue(res_s1_start)
            ca2 = get_ca_from_residue(res_s1_end)
            ca3 = get_ca_from_residue(res_s2_start)
            ca4 = get_ca_from_residue(res_s2_end)
            cb_pre = get_cb_from_residue(res_s1_end)

            if None in (ca1, ca2, ca3, ca4, cb_pre):
                continue

            # Strand geometry vectors
            first_strand_axis = ca2 - ca1
            midpoint_strand1 = (ca1 + ca2) / 2.0
            midpoint_strand2 = (ca3 + ca4) / 2.0
            sidechain_vector = cb_pre - ca2

            handedness, handedness_magnitude = compute_handedness_from_strand_vectors(
                first_strand_axis, midpoint_strand1, midpoint_strand2, sidechain_vector)

            hairpin_sequence = ''.join(chain_df.loc[s1_start : s2_end, "AA"].tolist())
            loop_sequence = ''.join(chain_df.loc[s1_end + 1 : s2_start - 1, "AA"].tolist()) if s2_start > s1_end + 1 else ""

            hairpins.append({
                "PDB": pdb_id,
                "Chain": chain_id,
                "strand1_start": chain_df.loc[s1_start, "ResNum"],
                "strand1_end":   chain_df.loc[s1_end,   "ResNum"],
                "strand2_start": chain_df.loc[s2_start, "ResNum"],
                "strand2_end":   chain_df.loc[s2_end,   "ResNum"],
                "loop_len":      s2_start - s1_end - 1,
                "handedness":    handedness,
                "handedness_magnitude": handedness_magnitude,
                "FullChainSequence": chain_sequence,
                "HairpinSequence": hairpin_sequence,
                "LoopSequence": loop_sequence
            })

    return pd.DataFrame(hairpins)


def compute_handedness_from_backbone_geometry(residue_pre, residue_post):
    """
    Compute ββ-hairpin handedness using backbone atom vectors as in Methods Summary.

    - u: N → C of residue before loop
    - v: CA_pre → CA_post
    - n: CA → CB of residue before loop

    Returns:
        handedness ("L"/"R") and normalized scalar triple product
    """
    if not all(atom in residue_pre for atom in ("N", "CA", "C")) or "CA" not in residue_post:
        return None, None

    u = Vector(residue_pre["C"].get_coord()) - Vector(residue_pre["N"].get_coord())
    v = Vector(residue_post["CA"].get_coord()) - Vector(residue_pre["CA"].get_coord())

    cb = get_cb_from_residue(residue_pre)
    ca = get_ca_from_residue(residue_pre)
    if cb is None or ca is None:
        return None, None

    n = cb - ca
    u_arr, v_arr, n_arr = u.get_array(), v.get_array(), n.get_array()
    numerator = np.dot(np.cross(u_arr, v_arr), n_arr)
    denominator = np.linalg.norm(u_arr) * np.linalg.norm(v_arr) * np.linalg.norm(n_arr)
    magnitude = numerator / denominator if denominator > 1e-6 else 0.0
    handedness = "R" if magnitude > 0 else "L"

    return handedness, magnitude


def detect_hairpins_and_chirality_backbone(model, dssp_df):
    """
    Detect β-hairpins and compute handedness using backbone-based vector definition.

    This follows the description from the Methods Summary:
    - Chirality is determined using N→C and CA→CA vectors around the loop
    - More robust to twisted β-strands

    Parameters:
        model    : Bio.PDB Model object
        dssp_df  : DSSP output DataFrame

    Returns:
        pd.DataFrame with per-hairpin handedness and sequence metadata
    """
    pdb_id = getattr(getattr(model, "parent", None), "id", "UNKNOWN")
    hairpins = []

    for chain_id in dssp_df["Chain"].unique():
        chain_df = dssp_df[dssp_df["Chain"] == chain_id].reset_index(drop=True)
        strands = identify_strands(chain_df)
        chain_sequence = ''.join(chain_df["AA"].tolist())

        for i in range(len(strands) - 1):
            s1_start, s1_end = strands[i]
            s2_start, s2_end = strands[i + 1]
            if np.any(chain_df.loc[s1_end + 1 : s2_start - 1, "SimpleSS"] != "C"):
                continue

            loop_len = s2_start - s1_end - 1
            if loop_len < 0:
                continue

            try:
                res_pre = model[chain_id][(" ", chain_df.loc[s1_end, "ResNum"], " ")]
                res_post = model[chain_id][(" ", chain_df.loc[s2_start, "ResNum"], " ")]
            except KeyError:
                continue

            handedness, magnitude = compute_handedness_from_backbone_geometry(res_pre, res_post)
            if handedness is None:
                continue

            hairpin_sequence = ''.join(chain_df.loc[s1_start : s2_end, "AA"].tolist())
            loop_sequence = ''.join(chain_df.loc[s1_end + 1 : s2_start - 1, "AA"].tolist()) if loop_len > 0 else ""

            hairpins.append({
                "PDB": pdb_id,
                "Chain": chain_id,
                "strand1_start": chain_df.loc[s1_start, "ResNum"],
                "strand1_end":   chain_df.loc[s1_end,   "ResNum"],
                "strand2_start": chain_df.loc[s2_start, "ResNum"],
                "strand2_end":   chain_df.loc[s2_end,   "ResNum"],
                "loop_len":      loop_len,
                "handedness":    handedness,
                "handedness_magnitude": magnitude,
                "FullChainSequence": chain_sequence,
                "HairpinSequence": hairpin_sequence,
                "LoopSequence": loop_sequence
            })

    return pd.DataFrame(hairpins)
