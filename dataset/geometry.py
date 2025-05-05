from Bio.PDB import Vector
import numpy as np


def get_ca_from_residue(residue):
    """Return the Cα Vector from a residue or None if missing."""
    if "CA" in residue:
        return Vector(residue["CA"].get_coord())
    return None




def get_cb_from_residue(residue):
    """
    Return Cβ as a Bio.PDB.Vector.

    • If a real “CB” atom exists → return it.
    • Else build an idealised Cβ from backbone N–CA–C:

          n̂ = unit(CA→N)                109.5°
          ĉ = unit(CA→C)             N\         /C
          n⊥ = n̂ − (n̂·ĉ)ĉ               \     /
          dir = –½ n⊥ + √3⁄2 (ĉ×n⊥)          CA ●───► Cβ*
          Cβ* = CA + 1.522 Å · dir

      (dir is 109.5 ° from both backbone bonds, giving a true tetrahedral
      side‑chain position.)

    Returns
    -------
    Bio.PDB.Vector
        Real or pseudo‑Cβ coordinates, or *None* if backbone atoms are missing.
    """
    # ───────────────────────────────────────────────────────────────── real CB
    if "CB" in residue:
        return Vector(residue["CB"].get_coord())

    # ───────────────────────────────────────────────────── require backbone N,C
    if not all(a in residue for a in ("N", "CA", "C")):
        return None

    N  = Vector(residue["N"].get_coord())
    CA = Vector(residue["CA"].get_coord())
    C  = Vector(residue["C"].get_coord())

    # numpy arrays for safe scalar maths
    n_hat = (N - CA).normalized().get_array()
    c_hat = (C - CA).normalized().get_array()

    # component of n̂ orthogonal to ĉ
    n_perp = n_hat - (n_hat @ c_hat) * c_hat
    n_perp /= np.linalg.norm(n_perp)

    # +120° rotation about ĉ :  dir = –½ n_perp + √3/2 (ĉ × n_perp)
    dir = -0.5 * n_perp + (np.sqrt(3) / 2.0) * np.cross(c_hat, n_perp)
    dir /= np.linalg.norm(dir)

    pseudo_cb = CA.get_array() + CA_CB_BOND * dir
    return Vector(pseudo_cb)


def identify_strands(df, min_len=2):
    segments = []
    in_seg = False
    for i, row in df.iterrows():
        if row["SimpleSS"] == "E":
            if not in_seg:
                start = i
                in_seg = True
        else:
            if in_seg and i - start >= min_len:
                segments.append((start, i - 1))
            in_seg = False
    if in_seg and len(df) - start >= min_len:
        segments.append((start, len(df) - 1))
    return segments
