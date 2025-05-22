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


## 3 · Dataset structure

| idx | PDB  | Chain | strand1_start | strand1_end | strand2_start | strand2_end | loop_len | handedness | handedness_mag | FullChainSequence          | HairpinSequence                   | LoopSequence   |
|----:|:----:|:-----:|--------------:|------------:|--------------:|------------:|---------:|:----------:|---------------:|:---------------------------|:----------------------------------|:--------------|
| 0 | 3vor | A | 69  | 70  | 75  | 76  | 4 | **L** | −0.555 | GSDSRTVSE… | RNGISGDY | GISG |
| 1 | 3vor | A | 75  | 76  | 78  | 83  | 1 | **L** | −0.622 | GSDSRTVSE… | DYIGIGGAI | *I* |
| 2 | 3vor | A | 78  | 83  | 93  | 101 | 9 | **L** | −0.681 | GSDSRTVSE… | GIGGAITSSGSTINKGFAMELNGL | TSSGSTINK |
| 3 | 3vor | A | 138 | 139 | 149 | 151 | 9 | **L** | −0.786 | GSDSRTVSE… | VNMLAATDNTTILR | MLAATDNTT |
| 4 | 5gji | A | 354 | 359 | 368 | 374 | 8 | **L** | −0.229 | SNMSLQNAE… | TFLVRDASTKMHGDYTLTLRK | ASTKMHGD |
| 5 | 5gji | A | 368 | 374 | 377 | 386 | 2 | **L** | −0.244 | SNMSLQNAE… | YTLTLRKGGNNKLIKIFHR | GG |
| 6 | 5gji | A | 377 | 386 | 389 | 391 | 2 | **L** | −0.420 | SNMSLQNAE… | NNKLIKIFHRDGKYG | DG |
| 7 | 5gji | A | 389 | 391 | 398 | 399 | 6 | **L** | −0.275 | SNMSLQNAE… | KYGFSDPLTFS | FSDPLT |
| 8 | 4nsv | A | 8   | 9   | 25  | 30  | 15| **R** |  0.415 | GVSGSCNID… | IDVVCPEGNGHRDVIRSVAAYSR | VVCPEGNGHRDVIRS |
| 9 | 4nsv | A | 25  | 30  | 33  | 41  | 2 | **R** |  0.956 | GVSGSCNID… | VAAYSRQGTMWCTGSLV | QG |



## 4 · Plot of loop length versus handedness for all 11,116 high-quality pdb entries

Local method of loop chirality definition was used.

Additional filter used on abs(handedness_magnitude) > 0.75

<img width="344" alt="image" src="https://github.com/user-attachments/assets/ed70bc1e-0f1f-4cd6-a0eb-a134f26b644d" />





## Appendix · Exact vs Backbone ββ‑chirality methods OUTDATED!

| Function (in `dataset/motif_logic.py`) | Vector **u** | Vector **v** | Vector **n** | Best for |
|----------------------------------------|--------------|--------------|--------------|----------|
| `detect_hairpins_and_chirality` <br>*(“exact / strand‑axis”)* | First‑strand axis  (Cα<sub>end</sub> − Cα<sub>start</sub>) | Midpoint<sub>strand1</sub> → Midpoint<sub>strand2</sub> | Cα→Cβ of residue before the loop | Fast; matches most published β‑hairpin surveys | 
| `detect_hairpins_and_chirality_backbone` <br>*(“backbone”)* | Backbone **N→C** of residue before loop | Cα<sub>pre</sub> → Cα<sub>post</sub> (across the loop) | Same Cα→Cβ side‑chain vector | Robust near the loop; insensitive to global sheet twist | 

Both return two chirality columns:

* `handedness` → **“L”** or **“R”** (sign of the scalar triple product)  
* `handedness_magnitude` → |scalar triple product|, *normalised*, ∈ [-1 … 1]

A magnitude close to 0 indicates an almost planar ββ‑unit where handedness is
ambiguous.  Perhaps we should keep only “strong” events:


