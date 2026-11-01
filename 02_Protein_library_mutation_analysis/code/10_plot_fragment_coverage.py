#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.13"
# dependencies = ["numpy==2.2.6", "polars==2.0.0", "matplotlib==3.10.7"]
# ///
# ─── How to run ───
# Install uv: https://docs.astral.sh/uv/getting-started/installation/
# Run the analysis entry point: python -m uv run code/09_calculate_fragment_coverage.py
# ──────────────────
"""Export clearly labelled coverage curves and discrete detection heatmaps."""

from __future__ import annotations

from pathlib import Path
from typing import Final

import matplotlib
import numpy as np
import polars as pl

matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch

ORDER: Final = "RKHDESTNQCGPAVILMFYWX"


def draw_coverage_by_position(root: Path, positions: pl.DataFrame) -> None:
    """Export the Count1 coverage curve with Arial text and a clean Protein position axis."""
    output = root / "figures/library_coverage"
    output.mkdir(parents=True, exist_ok=True)
    x = np.asarray(positions["FASTA_Position"].to_list(), dtype=np.int64)
    available = np.asarray(positions["Has_Input_Records"].to_list(), dtype=np.bool_)
    y = np.asarray(positions["Count1_Standard_AA_Count"].to_list(), dtype=np.float64)
    with plt.rc_context({
        "font.family": "Arial",
        "font.size": 10,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.unicode_minus": False,
    }):
        fig, axis = plt.subplots(figsize=(15, 4.5), layout="constrained")
        axis.plot(x, np.where(available, y, np.nan), color="#24778a", lw=1.4)
        for missing in x[~available]:
            axis.axvspan(missing - 0.5, missing + 0.5, color="0.8", alpha=0.8)
        axis.set(
            ylabel="Detected amino acid types (n)",
            xlabel="Protein position",
            ylim=(0, 21),
            xlim=(x[0] - 0.5, x[-1] + 0.5),
        )
        axis.set_yticks(np.arange(0, 21, 5))
        axis.set_xticks([int(x[0]), *range(50, int(x[-1]), 50), int(x[-1])])
        axis.grid(axis="y", color="0.9", lw=0.6)
        axis.spines[["top", "right"]].set_visible(False)
        axis.text(
            0.025, 0.15,
            f"Mean detected: {y.mean():.2f} amino acid types\nCoverage: {y.sum() / (20 * len(y)):.2%}",
            transform=axis.transAxes, fontsize=12, color="#245363", linespacing=1.7,
        )
        for extension in ("png", "pdf"):
            fig.savefig(output / f"coverage_by_position.{extension}", dpi=300)
        plt.close(fig)


def draw_figures(root: Path, positions: pl.DataFrame, grid: pl.DataFrame) -> None:
    """Export the Count1 coverage curve and threshold-specific detection heatmaps."""
    output = root / "figures/library_coverage"
    output.mkdir(parents=True, exist_ok=True)
    draw_coverage_by_position(root, positions)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "pdf.fonttype": 42})
    x = np.asarray(positions["FASTA_Position"].to_list(), dtype=np.int64)
    available = np.asarray(positions["Has_Input_Records"].to_list(), dtype=np.bool_)

    palette = ["#ffffff", "#39889a", "#333333", "#c97132", "#d6d6d6"]
    colors = ListedColormap(palette)
    legend = [Patch(facecolor=color, edgecolor="0.7", label=label) for color, label in zip(
        palette, ("Not detected at threshold", "Standard AA detected", "WT detected", "Stop X detected", "Position has no input"), strict=True,
    )]
    for threshold in (1, 5):
        matrix = np.zeros((len(ORDER), len(x)), dtype=np.int64)
        matrix[:, ~available] = 4
        for row_index, aa in enumerate(ORDER):
            part = grid.filter(pl.col("Observed_AA") == aa).sort("FASTA_Position")
            detected = np.asarray(part["Max_Fragment_Count"].to_list(), dtype=np.int64) >= threshold
            wt = np.asarray(part["Is_Wildtype"].to_list(), dtype=np.bool_)
            matrix[row_index, detected] = 3 if aa == "X" else 1
            matrix[row_index, detected & wt] = 2
        fig, axes = plt.subplots(2, 1, figsize=(18, 8), layout="constrained")
        split = (len(x) + 1) // 2
        for axis, start, stop in zip(axes, (0, split), (split, len(x)), strict=True):
            first, last = int(x[start]), int(x[stop - 1])
            axis.imshow(matrix[:, start:stop], aspect="auto", interpolation="none", cmap=colors,
                        norm=BoundaryNorm(np.arange(-0.5, 5, 1), colors.N), extent=(first - 0.5, last + 0.5, len(ORDER) - 0.5, -0.5))
            axis.set_yticks(np.arange(len(ORDER)), labels=list(ORDER))
            axis.set_ylabel("Observed amino acid (X = stop)")
            ticks = sorted({first, last, *range(((first + 24) // 25) * 25, last + 1, 25)})
            axis.set_xticks(ticks)
            axis.set_xlabel("FASTA position (data position = FASTA - 7)")
            axis.axhline(len(ORDER) - 1.5, color="0.45", lw=0.7)
        fig.suptitle(f"Group-0 amino-acid detection | Count >= {threshold} in at least one fragment | FASTA 1-8 excluded", fontsize=13)
        fig.legend(handles=legend, loc="outside lower center", ncol=5, frameon=False, fontsize=9)
        for extension in ("png", "pdf"):
            fig.savefig(output / f"amino_acid_presence_count{threshold}.{extension}", dpi=300)
        plt.close(fig)
