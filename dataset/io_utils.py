import os
import requests
import warnings
import pandas as pd
from Bio import PDB
from Bio.PDB import MMCIFParser, DSSP
from dataset.constants import simplified_ss_map
from pathlib import Path

def fetch_mmcif_file(pdb_id: str, outdir: str = ".", overwrite: bool = False) -> str:
    """
    Download an mmCIF file for *pdb_id* only if it is not already on disk.

    Parameters
    ----------
    pdb_id : str
        4-letter PDB code (extra characters are ignored, case‐insensitive).
    outdir : str
        Directory where the file will be saved (created if missing).
    overwrite : bool, default False
        Re-download even if the file exists.

    Returns
    -------
    str
        Absolute path to the local mmCIF file.
    """
    pdb_id = pdb_id.lower()[:4]
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    path = outdir / f"{pdb_id}.cif"

    if path.exists() and not overwrite:
        return str(path)

    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    r = requests.get(url, timeout=30)
    r.raise_for_status()
    path.write_text(r.text)

    return str(path)

def run_dssp_on_mmcif(pdb_path: str, dssp_exe="mkdssp"):
    """
    Parse an mmCIF file and run DSSP.
    On DSSP failure, issue a warning and return (None, None).
    """
    pdb_id = os.path.basename(pdb_path).split(".")[0]
    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure(pdb_id, pdb_path)
    model = structure[0]

    try:
        dssp = DSSP(model, pdb_path, dssp=dssp_exe)
    except Exception as e:
        warnings.warn(f"DSSP failed for {pdb_id}: {e}", RuntimeWarning)
        return model, None   # ← soft‑fail

    rows = []
    for key in dssp.keys():
        chain_id, (hetatm_flag, resseq, icode) = key
        aa = dssp[key][1]
        sec_str = dssp[key][2]
        simp = simplified_ss_map.get(sec_str, "C")
        asa = dssp[key][3]
        phi = dssp[key][4]
        psi = dssp[key][5]
        rows.append((chain_id, resseq, aa, sec_str, simp, asa, phi, psi))

    cols = ["Chain", "ResNum", "AA", "SecStruct", "SimpleSS", "ASA", "Phi", "Psi"]
    return model, pd.DataFrame(rows, columns=cols)