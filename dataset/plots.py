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


def plot_loop_length_chirality(
    df: pd.DataFrame,
    loop_lengths: Sequence[int] = (2, 3, 4, 5),
    ax: Optional[plt.Axes] = None,
):
    """
    Bar‑plot the count of β‑hairpins by loop length and chirality.

    Parameters
    ----------
    df : pd.DataFrame
        Must contain columns "loop_len" and "handedness".
    loop_lengths : iterable of int, default (2,3,4,5)
        Loop lengths to include on the x‑axis (order preserved).
    ax : matplotlib.axes.Axes, optional
        Draw into this axis; if None, create a new figure+axis.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4))

    # group counts
    subset = df[df["loop_len"].isin(loop_lengths)]
    counts = (
        subset.groupby(["loop_len", "handedness"])
               .size()
               .unstack(fill_value=0)
               .reindex(loop_lengths)
    )

    bar_w = 0.35
    x = range(len(counts))

    # L – black
    ax.bar([i - bar_w/2 for i in x],
           counts.get("L", 0),
           bar_w,
           label="L",
           edgecolor="black",
           facecolor="black")
    # R – white with black edge
    ax.bar([i + bar_w/2 for i in x],
           counts.get("R", 0),
           bar_w,
           label="R",
           edgecolor="black",
           facecolor="white")

    ax.set_xlabel("Loop length")
    ax.set_ylabel("Count")
    ax.set_xticks(x)
    ax.set_xticklabels(loop_lengths)
    ax.set_title("β‑hairpin loop length vs chirality")
    ax.legend(frameon=False)
    ax.figure.tight_layout()
    return ax
