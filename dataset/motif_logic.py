import numpy as np
import pandas as pd
from Bio.PDB import Vector
from .constants import CA_CB_BOND
from .geometry import (
    get_ca_from_residue,
    get_cb_from_residue,
    evaluate_triple_product_handedness,
)

# ──────────────────────────────────────────────────────────────────────
#  SECONDARY-STRUCTURE SEGMENT FINDERS
# ──────────────────────────────────────────────────────────────────────

def identify_strands(dssp_df, min_len: int = 2):
    """Return list of (start_idx, end_idx) for β-strands."""
    strands, in_seg = [], False
    for i, row in dssp_df.iterrows():
        if row["SimpleSS"] == "E":                       # β-strand
            if not in_seg:
                seg_start, in_seg = i, True
        else:                                            # exit strand
            if in_seg and (i - seg_start) >= min_len:
                strands.append((seg_start, i - 1))
            in_seg = False
    if in_seg and (len(dssp_df) - seg_start) >= min_len:
        strands.append((seg_start, len(dssp_df) - 1))
    return strands


def identify_helices(dssp_df, min_len: int = 5):
    """Return list of (start_idx, end_idx) for α-helices."""
    helices, in_seg = [], False
    for i, row in dssp_df.iterrows():
        if row["SimpleSS"] == "H":                       # α-helix
            if not in_seg:
                seg_start, in_seg = i, True
        else:
            if in_seg and (i - seg_start) >= min_len:
                helices.append((seg_start, i - 1))
            in_seg = False
    if in_seg and (len(dssp_df) - seg_start) >= min_len:
        helices.append((seg_start, len(dssp_df) - 1))
    return helices


# ──────────────────────────────────────────────────────────────────────
#  MOTIF DETECTION
# ──────────────────────────────────────────────────────────────────────

def detect_hairpins(dssp_df):
    """ββ hairpins with loop ≤ 5 residues."""
    hairpins = []
    for chain_id in dssp_df["Chain"].unique():
        chain_df = dssp_df[dssp_df["Chain"] == chain_id].reset_index(drop=True)
        strands = identify_strands(chain_df, min_len=2)
        chain_sequence = "".join(chain_df["AA"].tolist())

        for i in range(len(strands) - 1):
            s1_start, s1_end = strands[i]
            s2_start, s2_end = strands[i + 1]

            # ensure coil-only loop
            if np.any(chain_df.loc[s1_end + 1 : s2_start - 1, "SimpleSS"] != "C"):
                continue
            loop_len = s2_start - s1_end - 1
            if loop_len < 0 or loop_len > 5:
                continue

            hairpin_sequence = "".join(chain_df.loc[s1_start : s2_end, "AA"].tolist())
            loop_sequence = (
                "".join(chain_df.loc[s1_end + 1 : s2_start - 1, "AA"].tolist())
                if loop_len > 0
                else ""
            )

            hairpins.append(
                {
                    "Chain": chain_id,
                    "strand1_start_idx": s1_start,
                    "strand1_end_idx": s1_end,
                    "strand2_start_idx": s2_start,
                    "strand2_end_idx": s2_end,
                    "strand1_start_res": chain_df.loc[s1_start, "ResNum"],
                    "strand1_end_res": chain_df.loc[s1_end, "ResNum"],
                    "strand2_start_res": chain_df.loc[s2_start, "ResNum"],
                    "strand2_end_res": chain_df.loc[s2_end, "ResNum"],
                    "loop_len": loop_len,
                    "HairpinSequence": hairpin_sequence,
                    "LoopSequence": loop_sequence,
                    "FullChainSequence": chain_sequence,
                }
            )
    return hairpins


def _collect_backbone_vector(residues, n_atoms: int = 11):
    """Average coordinate of up to *n_atoms* backbone atoms (N, C, CA)."""
    coords = []
    for res in residues:
        for atom_name in ("N", "CA", "C"):
            if atom_name in res:
                coords.append(res[atom_name].get_coord())
            if len(coords) >= n_atoms:
                break
        if len(coords) >= n_atoms:
            break
    return Vector(np.mean(coords, axis=0)) if coords else None


def detect_ab_motifs(dssp_df, max_loop_len: int = 5):
    """αβ units (helix → strand)."""
    motifs = []
    for chain_id in dssp_df["Chain"].unique():
        chain_df = dssp_df[dssp_df["Chain"] == chain_id].reset_index(drop=True)
        helices = identify_helices(chain_df, min_len=3)
        strands = identify_strands(chain_df, min_len=2)
        chain_sequence = "".join(chain_df["AA"].tolist())

        for h_start, h_end in helices:
            after = [s for s in strands if s[0] > h_end]
            if not after:
                continue
            s_start, s_end = after[0]

            if np.any(chain_df.loc[h_end + 1 : s_start - 1, "SimpleSS"] != "C"):
                continue
            loop_len = s_start - h_end - 1
            if loop_len < 0 or loop_len > max_loop_len:
                continue

            motifs.append(
                {
                    "Chain": chain_id,
                    "helix_start_idx": h_start,
                    "helix_end_idx": h_end,
                    "strand_start_idx": s_start,
                    "strand_end_idx": s_end,
                    "helix_start_res": chain_df.loc[h_start, "ResNum"],
                    "helix_end_res": chain_df.loc[h_end, "ResNum"],
                    "strand_start_res": chain_df.loc[s_start, "ResNum"],
                    "strand_end_res": chain_df.loc[s_end, "ResNum"],
                    "loop_len": loop_len,
                    "FullChainSequence": chain_sequence,
                }
            )
    return motifs


def detect_ba_motifs(dssp_df, max_loop_len: int = 5):
    """βα units (strand → helix)."""
    motifs = []
    for chain_id in dssp_df["Chain"].unique():
        chain_df = dssp_df[dssp_df["Chain"] == chain_id].reset_index(drop=True)
        strands = identify_strands(chain_df, min_len=2)
        helices = identify_helices(chain_df, min_len=3)
        chain_sequence = "".join(chain_df["AA"].tolist())

        for s_start, s_end in strands:
            after = [h for h in helices if h[0] > s_end]
            if not after:
                continue
            h_start, h_end = after[0]

            if np.any(chain_df.loc[s_end + 1 : h_start - 1, "SimpleSS"] != "C"):
                continue
            loop_len = h_start - s_end - 1
            if loop_len < 0 or loop_len > max_loop_len:
                continue

            motifs.append(
                {
                    "Chain": chain_id,
                    "strand_start_idx": s_start,
                    "strand_end_idx": s_end,
                    "helix_start_idx": h_start,
                    "helix_end_idx": h_end,
                    "strand_start_res": chain_df.loc[s_start, "ResNum"],
                    "strand_end_res": chain_df.loc[s_end, "ResNum"],
                    "helix_start_res": chain_df.loc[h_start, "ResNum"],
                    "helix_end_res": chain_df.loc[h_end, "ResNum"],
                    "loop_len": loop_len,
                    "FullChainSequence": chain_sequence,
                }
            )
    return motifs


# ──────────────────────────────────────────────────────────────────────
#  ββ CHIRALITY ASSIGNERS
# ──────────────────────────────────────────────────────────────────────

def assign_beta_chirality_strand_axis(model, dssp_df, hairpin_annotations):
    """Triple-product handedness using strand-axis definition (ββ)."""
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
            continue  # missing residues in structure

        ca1, ca2 = get_ca_from_residue(res_s1_start), get_ca_from_residue(res_s1_end)
        ca3, ca4 = get_ca_from_residue(res_s2_start), get_ca_from_residue(res_s2_end)
        cb_pre   = get_cb_from_residue(res_s1_end)

        if None in (ca1, ca2, ca3, ca4, cb_pre):
            continue

        u = ca2 - ca1
        v = (ca4 + ca3) / 2.0 - (ca2 + ca1) / 2.0
        n = cb_pre - ca2
        handedness, magnitude = evaluate_triple_product_handedness(u, v, n)

        # copy everything except internal *_idx bookkeeping
        base = {k: v for k, v in h.items() if not k.endswith("_idx")}
        base.update(
            {
                "PDB": pdb_id,
                "handedness": handedness,
                "handedness_magnitude": magnitude,
            }
        )
        results.append(base)

    return pd.DataFrame(results)


def assign_beta_chirality_local(model, dssp_df, hairpin_annotations):
    """Triple-product handedness using local-geometry definition (ββ)."""
    pdb_id, results = getattr(getattr(model, "parent", None), "id", "UNKNOWN"), []

    for h in hairpin_annotations:
        chain_id = h["Chain"]
        chain_df = dssp_df[dssp_df["Chain"] == chain_id].reset_index(drop=True)

        try:
            res_pre  = model[chain_id][(" ", chain_df.loc[h["strand1_end_idx"], "ResNum"], " ")]
            res_post = model[chain_id][(" ", chain_df.loc[h["strand2_start_idx"], "ResNum"], " ")]
        except KeyError:
            continue

        if not all(atom in res_pre for atom in ("N", "CA", "C")) or "CA" not in res_post:
            continue

        cb, ca = get_cb_from_residue(res_pre), get_ca_from_residue(res_pre)
        if cb is None or ca is None:
            continue

        u = Vector(res_pre["C"].get_coord())  - Vector(res_pre["N"].get_coord())
        v = Vector(res_post["CA"].get_coord()) - Vector(res_pre["CA"].get_coord())
        n = cb - ca
        handedness, magnitude = evaluate_triple_product_handedness(u, v, n)

        base = {k: v for k, v in h.items() if not k.endswith("_idx")}
        base.update(
            {
                "PDB": pdb_id,
                "handedness": handedness,
                "handedness_magnitude": magnitude,
            }
        )
        results.append(base)

    return pd.DataFrame(results)


# ──────────────────────────────────────────────────────────────────────
#  αβ / βα ORIENTATION ASSIGNERS
# ──────────────────────────────────────────────────────────────────────

def _orientation_class(u: Vector, n: Vector, eps: float = 1e-6):
    """Return ('P'|'A', |cosθ|) where θ is angle between *u* and *n*."""
    if u.norm() < eps or n.norm() < eps:
        return None, None
    angle_rad = u.angle(n)
    orientation = "P" if angle_rad < (np.pi / 2) else "A"
    return orientation, abs(np.cos(angle_rad))


def assign_ab_chirality(model, dssp_df, ab_annotations):
    """Orientation of αβ units (helix → strand)."""
    pdb_id, results = getattr(getattr(model, "parent", None), "id", "UNKNOWN"), []

    for m in ab_annotations:
        chain_id = m["Chain"]
        chain_df = dssp_df[dssp_df["Chain"] == chain_id].reset_index(drop=True)

        try:
            helix_end_res = model[chain_id][(" ", chain_df.loc[m["helix_end_idx"], "ResNum"], " ")]
            strand_start_res = model[chain_id][(" ", chain_df.loc[m["strand_start_idx"], "ResNum"], " ")]
        except KeyError:
            continue

        # avg of last 11 backbone atoms in helix
        helix_residues = [
            model[chain_id][(" ", chain_df.loc[idx, "ResNum"], " ")]
            for idx in range(max(m["helix_end_idx"] - 3, m["helix_start_idx"]), m["helix_end_idx"] + 1)
        ]
        avg_helix_vec = _collect_backbone_vector(helix_residues)
        ca_strand, cb_strand = get_ca_from_residue(strand_start_res), get_cb_from_residue(strand_start_res)
        if None in (avg_helix_vec, ca_strand, cb_strand):
            continue

        u, n = ca_strand - avg_helix_vec, cb_strand - ca_strand
        orient, mag = _orientation_class(u, n)
        if orient is None:
            continue

        base = {k: v for k, v in m.items() if not k.endswith("_idx")}
        base.update({"PDB": pdb_id, "orientation": orient, "orientation_magnitude": mag})
        results.append(base)

    return pd.DataFrame(results)


def assign_ba_chirality(model, dssp_df, ba_annotations):
    """Orientation of βα units (strand → helix)."""
    pdb_id, results = getattr(getattr(model, "parent", None), "id", "UNKNOWN"), []

    for m in ba_annotations:
        chain_id = m["Chain"]
        chain_df = dssp_df[dssp_df["Chain"] == chain_id].reset_index(drop=True)

        try:
            strand_end_res = model[chain_id][(" ", chain_df.loc[m["strand_end_idx"], "ResNum"], " ")]
            helix_start_res = model[chain_id][(" ", chain_df.loc[m["helix_start_idx"], "ResNum"], " ")]
        except KeyError:
            continue

        # avg of first 11 backbone atoms in helix
        helix_residues = [
            model[chain_id][(" ", chain_df.loc[idx, "ResNum"], " ")]
            for idx in range(m["helix_start_idx"], min(m["helix_start_idx"] + 4, m["helix_end_idx"] + 1))
        ]
        avg_helix_vec = _collect_backbone_vector(helix_residues)
        ca_strand, cb_strand = get_ca_from_residue(strand_end_res), get_cb_from_residue(strand_end_res)
        if None in (avg_helix_vec, ca_strand, cb_strand):
            continue

        u, n = avg_helix_vec - ca_strand, cb_strand - ca_strand
        orient, mag = _orientation_class(u, n)
        if orient is None:
            continue

        base = {k: v for k, v in m.items() if not k.endswith("_idx")}
        base.update({"PDB": pdb_id, "orientation": orient, "orientation_magnitude": mag})
        results.append(base)

    return pd.DataFrame(results)
