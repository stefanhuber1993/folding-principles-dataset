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
    """Nature‑style bar chart of β‑hairpin chirality vs loop length."""
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots(figsize=(3, 2))

    # counts ---------------------------------------------------------------
    subset = df[df["loop_len"].isin(loop_lengths)]
    counts = (
        subset.groupby(["loop_len", "handedness"])
               .size()
               .unstack(fill_value=0)
               .reindex(loop_lengths)
    )

    bar_w = 0.35
    x = range(len(counts))

    # style ----------------------------------------------------------------
    ax.grid(axis="y", color="0.85", lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    # bars -----------------------------------------------------------------
    ax.bar([i - bar_w/2 for i in x],
           counts.get("L", 0),
           bar_w,
           color="black",
           edgecolor="black",
           linewidth=1.0,
           zorder=2)

    ax.bar([i + bar_w/2 for i in x],
           counts.get("R", 0),
           bar_w,
           color="white",
           edgecolor="black",
           linewidth=1.0,
           zorder=2)

    # text labels ----------------------------------------------------------
    # for i, length in enumerate(loop_lengths):
    #     for offset, hand in [(-bar_w/2, "L"), (bar_w/2, "R")]:
    #         val = counts.at[length, hand] if hand in counts.columns else 0
    #         if val == 0:
    #             continue
    #         ax.text(i + offset, val + 20, str(val),
    #                 ha="center", va="bottom",
    #                 fontsize=7)

    # axes -----------------------------------------------------------------
    ax.set_xticks(x)
    ax.set_xticklabels(loop_lengths, fontsize=9)
    ax.set_xlabel("Loop length", fontsize=9)
    ax.set_ylabel("Frequency", fontsize=9)
    ax.tick_params(axis="both", labelsize=8)
    ax.set_ylim(0, counts.values.max() * 1.15)

    return ax
