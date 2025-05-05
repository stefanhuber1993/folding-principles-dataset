from Bio.PDB import Vector

def get_ca_from_residue(residue):
    """Return the Cα Vector from a residue or None if missing."""
    if "CA" in residue:
        return Vector(residue["CA"].get_coord())
    return None


def get_cb_from_residue(residue):
    """
    Return a Vector for the Cβ atom of a residue.

    • If CB exists           → use it.
    • Else if backbone‑H     → use one of those.
    • Else                   → construct pseudo‑CB from N–CA–C geometry.

    Always returns a Vector or None.
    """
    # 1. real Cβ
    if "CB" in residue:
        return Vector(residue["CB"].get_coord())

    # 2. any H attached to CA
    for h in ("H", "HA", "1H", "2H", "HA2", "HA3"):
        if h in residue:
            return Vector(residue[h].get_coord())

    # 3. pseudo‑CB from backbone
    for atom in ("N", "CA", "C"):
        if atom not in residue:
            return None

    N  = Vector(residue["N"].get_coord())
    CA = Vector(residue["CA"].get_coord())
    C  = Vector(residue["C"].get_coord())

    v1 = (N - CA).normalized()
    v2 = (C - CA).normalized()

    # bisector (pointing away from backbone), then normalise
    direction = -(v1 + v2).normalized()            # still a Vector

    # scale the *array* explicitly → 1.522 Å
    direction_scaled = direction.get_array() * 1.522
    pseudo_cb = Vector(CA.get_array() + direction_scaled)

    return pseudo_cb


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
