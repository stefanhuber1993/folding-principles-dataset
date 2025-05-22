# dataset/plots.py
"""
Plot utilities for Folding-Principles Dataset
-------------------------------------------
Helper:

    plot_loop_length_category(df, loop_lengths=[2,3,4,5], ax=None,
                              category_col="handedness")

`df` must contain:
    • "loop_len"              – integer loop length
    • *category_col*          – one of
          - "handedness"  → category values "L", "R"
          - "orientation" → category values "P", "A"

Example
-------
>>> from dataset.plots import plot_loop_length_category
>>> plot_loop_length_category(all_hairpins)                # L / R
>>> plot_loop_length_category(all_ab_motifs, category_col="orientation")  # P / A
"""

from typing import Optional, Sequence, List
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# consistent order and simple black/white/grey style -------------------
_CATEGORY_ORDER: List[str] = ["L", "R", "P", "A"]
_FACECOLORS = {
    "L": "black",
    "R": "white",
    "P": "0.5",   # mid-grey
    "A": "white",
}
_EDGECOLORS = {
    "L": "black",
    "R": "black",
    "P": "black",
    "A": "black",
}
_HATCHES = {
    "L": None,
    "R": None,
    "P": None,
    "A": "//",    # hatched to distinguish from R
}


def plot_loop_length_category(
    df: pd.DataFrame,
    loop_lengths: Sequence[int] = (2, 3, 4, 5),
    *,
    category_col: str = "handedness",
    ax: Optional[plt.Axes] = None,
):
    """Nature-style bar chart of category frequencies vs β-hairpin loop length.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing at least "loop_len" and *category_col*.
    loop_lengths : sequence of int
        Loop lengths to display on the x-axis.
    category_col : {"handedness", "orientation"}
        Which column to use for the categories. The function adjusts to any
        combination of the letters L/R and P/A found in the data.
    ax : matplotlib.axes.Axes, optional
        Use an existing axes; otherwise a new one is created.

    Returns
    -------
    matplotlib.axes.Axes
    """
    if df.empty:
        print("Warning: empty DataFrame passed to plot_loop_length_category.")
        return ax

    if category_col not in df.columns:
        raise ValueError(f"Column '{category_col}' not found in DataFrame.")

    if ax is None:
        _, ax = plt.subplots(figsize=(2.2, 2.2), dpi=150)

    # ------------------------------------------------------------------
    # counts
    subset = df[df["loop_len"].isin(loop_lengths)]
    counts = (
        subset.groupby(["loop_len", category_col])
               .size()
               .unstack(fill_value=0)
               .reindex(loop_lengths, fill_value=0)
    )

    categories_present = [c for c in _CATEGORY_ORDER if c in counts.columns]
    n_cat = len(categories_present)
    if n_cat == 0:
        print("Warning: no recognised categories in the data.")
        return ax

    # ------------------------------------------------------------------
    # styling
    ax.grid(axis="y", color="0.85", lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)

    # ------------------------------------------------------------------
    # bars
    full_w = 0.8
    bar_w = full_w / n_cat
    x_pos = np.arange(len(loop_lengths))

    for i, cat in enumerate(categories_present):
        offs = (-full_w / 2) + (i + 0.5) * bar_w
        values = counts.get(cat, pd.Series([0] * len(loop_lengths)))
        bars = ax.bar(
            x_pos + offs,
            values,
            bar_w,
            facecolor=_FACECOLORS[cat],
            edgecolor=_EDGECOLORS[cat],
            linewidth=1.0,
            hatch=_HATCHES[cat],
            zorder=2,
            label=cat,
        )

    # ------------------------------------------------------------------
    # axes
    ax.set_xticks(x_pos)
    ax.set_xticklabels(loop_lengths, fontsize=9)
    ax.set_xlabel("Loop length", fontsize=9)
    ax.set_ylabel("Frequency", fontsize=9)
    ax.tick_params(axis="both", labelsize=8)

    # y-limit safeguard
    max_val = counts.values.max()
    ax.set_ylim(0, max_val * 1.15 if np.isfinite(max_val) and max_val > 0 else 1)

    # legend only for categories present
    ax.legend(frameon=False, fontsize=8, loc="upper right")

    return ax
