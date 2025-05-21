import numpy as np
import pandas as pd
from Bio.PDB import Vector
from .constants import CA_CB_BOND
from .geometry import get_ca_from_residue, get_cb_from_residue

def identify_strands(dssp_df, min_len=2):
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

def detect_hairpins(dssp_df):
    hairpins = []
    for chain_id in dssp_df["Chain"].unique():
        chain_df = dssp_df[dssp_df["Chain"] == chain_id].reset_index(drop=True)
        strands = identify_strands(chain_df, min_len=2)
        chain_sequence = ''.join(chain_df["AA"].tolist())

        for i in range(len(strands) - 1):
            s1_start, s1_end = strands[i]
            s2_start, s2_end = strands[i + 1]

            if np.any(chain_df.loc[s1_end + 1 : s2_start - 1, "SimpleSS"] != "C"):
                continue

            loop_len = s2_start - s1_end - 1
            if loop_len < 0:
                continue

            hairpin_sequence = ''.join(chain_df.loc[s1_start : s2_end, "AA"].tolist())
            loop_sequence = ''.join(chain_df.loc[s1_end + 1 : s2_start - 1, "AA"].tolist()) if loop_len > 0 else ""

            hairpins.append({
                "Chain": chain_id,
                "strand1_start_idx": s1_start,
                "strand1_end_idx": s1_end,
                "strand2_start_idx": s2_start,
                "strand2_end_idx": s2_end,
                "strand1_start_res": chain_df.loc[s1_start, "ResNum"],
                "strand1_end_res":   chain_df.loc[s1_end,   "ResNum"],
                "strand2_start_res": chain_df.loc[s2_start, "ResNum"],
                "strand2_end_res":   chain_df.loc[s2_end,   "ResNum"],
                "loop_len":          loop_len,
                "HairpinSequence":   hairpin_sequence,
                "LoopSequence":      loop_sequence,
                "FullChainSequence": chain_sequence
            })
    return hairpins

def evaluate_triple_product_handedness(u, v, n):
    u_arr = u.get_array() if isinstance(u, Vector) else u
    v_arr = v.get_array() if isinstance(v, Vector) else v
    n_arr = n.get_array() if isinstance(n, Vector) else n

    numerator = np.dot(np.cross(u_arr, v_arr), n_arr)
    denominator = np.linalg.norm(u_arr) * np.linalg.norm(v_arr) * np.linalg.norm(n_arr)
    magnitude = numerator / denominator if denominator > 1e-6 else 0.0
    handedness = "R" if magnitude > 0 else "L"
    return handedness, magnitude

def assign_beta_chirality_strand_axis(model, dssp_df, hairpin_annotations):
    pdb_id = getattr(getattr(model, "parent", None), "id", "UNKNOWN")
    results = []

    for h in hairpin_annotations:
        chain_id = h["Chain"]
        chain_df = dssp_df[dssp_df["Chain"] == chain_id].reset_index(drop=True)

        try:
            res_s1_start = model[chain_id][(" ", chain_df.loc[h["strand1_start_idx"], "ResNum"], " ")]
            res_s1_end   = model[chain_id][(" ", chain_df.loc[h["strand1_end_idx"],   "ResNum"], " ")]
            res_s2_start = model[chain_id][(" ", chain_df.loc[h["strand2_start_idx"], "ResNum"], " ")]
            res_s2_end   = model[chain_id][(" ", chain_df.loc[h["strand2_end_idx"],   "ResNum"], " ")]
        except KeyError:
            continue

        ca1 = get_ca_from_residue(res_s1_start)
        ca2 = get_ca_from_residue(res_s1_end)
        ca3 = get_ca_from_residue(res_s2_start)
        ca4 = get_ca_from_residue(res_s2_end)
        cb_pre = get_cb_from_residue(res_s1_end)

        if None in (ca1, ca2, ca3, ca4, cb_pre):
            continue

        u = ca2 - ca1
        v = (ca4 + ca3) / 2.0 - (ca2 + ca1) / 2.0
        n = cb_pre - ca2

        handedness, magnitude = evaluate_triple_product_handedness(u, v, n)

        results.append({
            **h,
            "PDB": pdb_id,
            "handedness": handedness,
            "handedness_magnitude": magnitude,
        })

    return pd.DataFrame(results)

def assign_beta_chirality_local(model, dssp_df, hairpin_annotations):
    pdb_id = getattr(getattr(model, "parent", None), "id", "UNKNOWN")
    results = []

    for h in hairpin_annotations:
        chain_id = h["Chain"]
        chain_df = dssp_df[dssp_df["Chain"] == chain_id].reset_index(drop=True)

        try:
            res_pre = model[chain_id][(" ", chain_df.loc[h["strand1_end_idx"], "ResNum"], " ")]
            res_post = model[chain_id][(" ", chain_df.loc[h["strand2_start_idx"], "ResNum"], " ")]
        except KeyError:
            continue

        if not all(atom in res_pre for atom in ("N", "CA", "C")) or "CA" not in res_post:
            continue

        cb = get_cb_from_residue(res_pre)
        ca = get_ca_from_residue(res_pre)
        if cb is None or ca is None:
            continue

        u = Vector(res_pre["C"].get_coord()) - Vector(res_pre["N"].get_coord())
        v = Vector(res_post["CA"].get_coord()) - Vector(res_pre["CA"].get_coord())
        n = cb - ca

        handedness, magnitude = evaluate_triple_product_handedness(u, v, n)

        results.append({
            **h,
            "PDB": pdb_id,
            "handedness": handedness,
            "handedness_magnitude": magnitude,
        })

    return pd.DataFrame(results)