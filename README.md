# Folding Principles Dataset

This repository generates a curated dataset of protein β-hairpins and their handedness, based on structural data from the PDB. It uses DSSP for secondary structure assignment and computes chirality for each hairpin loop.

---

## 🧪 Environment Setup (Conda, recommended)

We use a Conda environment to manage all dependencies, including the DSSP binary.

### ✅ Recommended: Create Conda environment from `environment.yml`

```bash
# If you're on macOS with an M1/M2/M3 (ARM64) chip:
export CONDA_SUBDIR=osx-64

conda env create -f environment.yml
conda activate folding-dssp

# If you're on macOS with an M1/M2/M3 (ARM64) chip:
conda config --env --set subdir osx-64
```
