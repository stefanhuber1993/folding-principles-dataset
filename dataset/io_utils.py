import os
import requests
import warnings
import pandas as pd
from Bio import PDB
from Bio.PDB import MMCIFParser, DSSP
from dataset.constants import simplified_ss_map

def fetch_mmcif_file(pdb_id: str, outdir: str = ".", overwrite: bool = False) -> str:
    """
    Download mmCIF file for given PDB ID and return local path.
    """
    pdb_id = pdb_id.lower()[:4]
    filename = f"{pdb_id}.cif"
    path = os.path.join(outdir, filename)

    if os.path.exists(path) and not overwrite:
        return path

    url = f"https://files.rcsb.org/download/{filename}"
    r = requests.get(url, timeout=10)
    r.raise_for_status()

    with open(path, "w") as f:
        f.write(r.text)

    return path

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