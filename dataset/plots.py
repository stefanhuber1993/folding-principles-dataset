# dataset/plots.py
"""
Plot utilities for Folding‑Principles Dataset
-------------------------------------------
Currently contains one helper:

    plot_loop_length_chirality(df, loop_lengths=[2,3,4,5], ax=None)

`df` must have at least two columns:
    • "loop_len"      – integer loop length
    • "handedness"    – "L" or "R"

Example
-------
>>> from dataset.plots import plot_loop_length_chirality
>>> plot_loop_length_chirality(all_hairpins)
"""

from typing import Iterable, Optional, Sequence
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def plot_loop_length_chirality(
    df: pd.DataFrame,
    loop_lengths: Sequence[int] = (2, 3, 4, 5),
    ax: Optional[plt.Axes] = None,
):
    """Nature‑style bar chart of β‑hairpin chirality vs loop length."""
    if df.empty:
        print("Warning: empty DataFrame passed to plot_loop_length_chirality.")
        return ax

    if ax is None:
        fig, ax = plt.subplots(figsize=(2, 2), dpi=150)

    # counts ---------------------------------------------------------------
    subset = df[df["loop_len"].isin(loop_lengths)]
    counts = (
        subset.groupby(["loop_len", "handedness"])
               .size()
               .unstack(fill_value=0)
               .reindex(loop_lengths, fill_value=0)
    )

    bar_w = 0.35
    x = range(len(counts))

    # style ----------------------------------------------------------------
    ax.grid(axis="y", color="0.85", lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    # bars -----------------------------------------------------------------
    ax.bar([i - bar_w / 2 for i in x],
           counts.get("L", pd.Series([0]*len(loop_lengths))),
           bar_w,
           color="black",
           edgecolor="black",
           linewidth=1.0,
           zorder=2, label='L')

    ax.bar([i + bar_w / 2 for i in x],
           counts.get("R", pd.Series([0]*len(loop_lengths))),
           bar_w,
           color="white",
           edgecolor="black",
           linewidth=1.0,
           zorder=2, label="R")

    # axes -----------------------------------------------------------------
    ax.set_xticks(x)
    ax.set_xticklabels(loop_lengths, fontsize=9)
    ax.set_xlabel("Loop length", fontsize=9)
    ax.set_ylabel("Frequency", fontsize=9)
    ax.tick_params(axis="both", labelsize=8)

    # safeguard for ylim --------------------------------------------------
    max_val = counts.values.max()
    ax.set_ylim(0, max_val * 1.15 if np.isfinite(max_val) and max_val > 0 else 1)

    ax.legend(frameon=False, fontsize=8, loc="upper right")

    return ax
