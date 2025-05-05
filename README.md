# Folding Principles Dataset – ββ‑Hairpin Prototype (WIP)

Extracts β‑hairpin (ββ‑unit) motifs from PDB mmCIF structures and analyses the
relationship between **loop length** and **chirality** (L / R), reproducing the
trend first reported in **Koga _et al._ 2012, *Nature* 491:222‑227**.

Only the ββ‑pattern is implemented for now.  The remaining Koga motifs
(αα‑corner, βαβ‑unit, αβ‑motif, etc.) will be added in later commits.

---

## 1 · Environment (set‑up once)

```bash
git clone https://github.com/your‑name/FoldingPrinciplesDataset.git
cd FoldingPrinciplesDataset

# (x86 Linux / macOS Intel)
conda env create -f environment.yml
conda activate folding-dssp

# (macOS Apple Silicon) – the dssp wheel is x86 only
export CONDA_SUBDIR=osx-64
conda env create -f environment.yml
conda activate folding-dssp
conda config --env --set subdir osx-64
```


## 2 · Run the notebook

The current workflow is demonstrated in **`notebooks/test_pipeline_parts.ipynb`**.

1. **Open the notebook** after activating the Conda environment.  
2. **Edit the `pdb_ids` cell** (or leave it empty to let the PISCES fetcher
   choose high‑quality structures automatically).  
3. **Run all cells**. The notebook will  
   * download the required mmCIF files  
   * run DSSP on each structure  
   * detect β‑hairpins  
   * compute chirality with both methods  
   * produce the loop‑length × handedness bar chart


## 3 · Example of a preliminary plot using 400 pdb files

Additional filter used on abs(handedness_magnitude) > 0.8
<img width="329" alt="image" src="https://github.com/user-attachments/assets/55dd9c4c-0d84-4962-b08b-4820ea9e0803" />






## Appendix · Exact vs Backbone ββ‑chirality methods

| Function (in `dataset/motif_logic.py`) | Vector **u** | Vector **v** | Vector **n** | Best for | Caveats |
|----------------------------------------|--------------|--------------|--------------|----------|---------|
| `detect_hairpins_and_chirality` <br>*(“exact / strand‑axis”)* | First‑strand axis  (Cα<sub>end</sub> − Cα<sub>start</sub>) | Midpoint<sub>strand1</sub> → Midpoint<sub>strand2</sub> | Cα→Cβ of residue before the loop | Fast; matches most published β‑hairpin surveys | Becomes noisy when strands are strongly twisted |
| `detect_hairpins_and_chirality_backbone` <br>*(“backbone”)* | Backbone **N→C** of residue before loop | Cα<sub>pre</sub> → Cα<sub>post</sub> (across the loop) | Same Cα→Cβ side‑chain vector | Robust near the loop; insensitive to global sheet twist | Slightly slower (needs two residue look‑ups) |

Both return two chirality columns:

* `handedness` → **“L”** or **“R”** (sign of the scalar triple product)  
* `handedness_magnitude` → |scalar triple product|, *normalised*, ∈ [0 … 1]

### Recommended filter

A magnitude close to 0 indicates an almost planar ββ‑unit where handedness is
ambiguous.  Perhaps we should keep only “strong” events:

```python
strong = hairpins[hq_df["handedness_magnitude"].abs() > 0.75]
```


