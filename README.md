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

