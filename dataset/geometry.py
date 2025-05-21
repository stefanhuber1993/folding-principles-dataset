from Bio.PDB import Vector
import numpy as np
from dataset.constants import CA_CB_BOND 


def get_ca_from_residue(residue):
    """Return the Cα Vector from a residue or None if missing."""
    if "CA" in residue:
        return Vector(residue["CA"].get_coord())
    return None


def get_cb_from_residue(residue):
    """
    Return Cβ coordinates as Bio.PDB.Vector.

    • If a real “CB” atom exists → return it.  
    • Otherwise build an *ideal* Cβ that is 109.5 ° from both backbone bonds,
      and placed on the true side-chain side (opposite the carbonyl C).

    Construction (works for any ∠N-CA-C):
    ------------------------------------------------
      u = unit (CA→N)
      v = unit (CA→C)

      â = (u − v)/‖u − v‖
      b̂ = (u + v)/‖u + v‖
      ŵ = â × b̂

      dir = −⅓ b̂ + (2√2/3) ŵ
      if dir·v > 0 → flip dir (select side-chain lobe)

      Cβ* = CA + CA_CB_BOND · dir
    """
    # 1 ───────────────────────────  real CB present
    if "CB" in residue:
        return Vector(residue["CB"].get_coord())

    # 2 ───────────────────────────  backbone atoms required
    if not all(a in residue for a in ("N", "CA", "C")):
        return None

    N  = Vector(residue["N"].get_coord())
    CA = Vector(residue["CA"].get_coord())
    C  = Vector(residue["C"].get_coord())

    # unit backbone vectors (NumPy arrays)
    u = (N - CA).normalized().get_array()    # CA→N
    v = (C - CA).normalized().get_array()    # CA→C

    # orthonormal basis in N-CA-C plane
    a = u - v;  a /= np.linalg.norm(a)
    b = u + v;  b /= np.linalg.norm(b)
    w = np.cross(a, b);  w /= np.linalg.norm(w)

    # tetrahedral direction (109.47° from u and v)
    dir_vec = (-1/3) * b + (2.0 * np.sqrt(2) / 3.0) * w
    dir_vec /= np.linalg.norm(dir_vec)

    # ensure Cβ lies opposite the carbonyl C (true side-chain side)
    if np.dot(dir_vec, v) > 0:
        dir_vec = -dir_vec

    pseudo_cb = CA.get_array() + CA_CB_BOND * dir_vec
    return Vector(pseudo_cb)


def evaluate_triple_product_handedness(u, v, n):
    u_arr = u.get_array() if isinstance(u, Vector) else u
    v_arr = v.get_array() if isinstance(v, Vector) else v
    n_arr = n.get_array() if isinstance(n, Vector) else n

    numerator = np.dot(np.cross(u_arr, v_arr), n_arr)
    denominator = np.linalg.norm(u_arr) * np.linalg.norm(v_arr) * np.linalg.norm(n_arr)
    magnitude = numerator / denominator if denominator > 1e-6 else 0.0
    handedness = "R" if magnitude > 0 else "L"
    return handedness, magnitude


# def identify_strands(df, min_len=2):
#     segments = []
#     in_seg = False
#     for i, row in df.iterrows():
#         if row["SimpleSS"] == "E":
#             if not in_seg:
#                 start = i
#                 in_seg = True
#         else:
#             if in_seg and i - start >= min_len:
#                 segments.append((start, i - 1))
#             in_seg = False
#     if in_seg and len(df) - start >= min_len:
#         segments.append((start, len(df) - 1))
#     return segments
