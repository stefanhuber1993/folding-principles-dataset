from Bio.PDB import Vector

def get_ca(model, chain_id, resseq, icode=" "):
    try:
        return Vector(model[chain_id][(" ", resseq, icode)]["CA"].get_coord())
    except Exception:
        return None



def get_cb(model, chain_id, resseq, icode=" "):
    """
    Return a Bio.PDB.Vector for the Cβ atom of a residue.
    If Cβ is missing:
      1) Return any backbone-attached hydrogen (if available)
      2) Else approximate pseudo-Cβ from N–CA–C geometry.
    """
    try:
        residue = model[chain_id][(" ", resseq, icode)]
    except KeyError:
        return None

    if "CB" in residue:
        return Vector(residue["CB"].get_coord())

    for h in ("H", "HA", "1H", "2H", "HA2", "HA3"):
        if h in residue:
            return Vector(residue[h].get_coord())

    for atom in ("N", "CA", "C"):
        if atom not in residue:
            return None

    N = Vector(residue["N"].get_coord())
    CA = Vector(residue["CA"].get_coord())
    C = Vector(residue["C"].get_coord())

    v1 = (N - CA).normalized()
    v2 = (C - CA).normalized()
    direction = (v1 + v2).normalized().get_array()

    # Use numpy array math and reconstruct Vector
    pseudo = Vector(CA.get_array() - 1.522 * direction)
    return pseudo

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
