#!/usr/bin/env python3
"""Plot mutation frequencies, replicate agreement and amino-acid coverage."""
import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr, spearmanr

AA = list("ACDEFGHIKLMNPQRSTVWY")
KEYS = ["reference_aa_position", "wt_aa", "mut_aa"]


def load_frequency(path: Path, round_name: str) -> pd.DataFrame:
    data = pd.read_csv(path)
    required = set(KEYS + ["library_round", "count", "frequency"])
    if missing := required - set(data.columns):
        raise ValueError(f"Missing columns in {path}: {sorted(missing)}")
    data = data.loc[data.library_round == round_name].copy()
    if data.empty:
        raise ValueError(f"No rows for {round_name!r} in {path}")
    # Each amino acid may be represented by multiple synonymous codon rows.
    return data.groupby(KEYS, as_index=False, dropna=False)[["count", "frequency"]].sum()


def plot_frequency(data: pd.DataFrame, destination: Path) -> None:
    mutations = data.loc[(data.wt_aa != data.mut_aa) & data.mut_aa.isin(AA)]
    matrix = mutations.pivot_table(index="mut_aa", columns="reference_aa_position", values="frequency", aggfunc="sum").reindex(AA)
    if matrix.empty:
        raise ValueError("No standard amino-acid mutations available for heatmap")
    fig, ax = plt.subplots(figsize=(max(8, min(24, len(matrix.columns) / 12)), 5))
    im = ax.imshow(matrix.to_numpy(dtype=float), aspect="auto", interpolation="nearest", cmap="viridis")
    ax.set_yticks(range(len(AA)), AA); ax.set_xlabel("Reference amino-acid position"); ax.set_ylabel("Mutant amino acid")
    ticks = np.linspace(0, len(matrix.columns) - 1, min(10, len(matrix.columns)), dtype=int)
    ax.set_xticks(ticks, [str(matrix.columns[i]) for i in ticks]); fig.colorbar(im, ax=ax, label="Mutation frequency")
    fig.tight_layout(); fig.savefig(destination, dpi=200); plt.close(fig)


def plot_coverage(first: pd.DataFrame, second: pd.DataFrame, destination: Path, table: Path) -> None:
    frames = []
    for label, data in (("replicate_1", first), ("replicate_2", second)):
        data = data.loc[(data.wt_aa != data.mut_aa) & data.mut_aa.isin(AA)].copy()
        for threshold in (1, 5):
            counts = data.loc[data['count'] >= threshold].groupby('reference_aa_position').mut_aa.nunique()
            positions = sorted(data.reference_aa_position.unique())
            frames.extend({"replicate": label, "threshold": threshold, "position": pos, "detected_non_wt": int(counts.get(pos, 0))} for pos in positions)
    result = pd.DataFrame(frames); result.to_csv(table, index=False)
    fig, ax = plt.subplots(figsize=(10, 4))
    for (replicate, threshold), group in result.groupby(["replicate", "threshold"]):
        ax.plot(group.position, group.detected_non_wt, label=f"{replicate}, count ≥ {threshold}")
    ax.set(xlabel="Reference amino-acid position", ylabel="Detected non-WT amino acids", ylim=(0, 19)); ax.legend(); fig.tight_layout(); fig.savefig(destination, dpi=200); plt.close(fig)


def enrichment_from_counts(path: Path, comparison: str, pseudocount: float) -> pd.DataFrame:
    later, separator, earlier = comparison.partition("_vs_")
    if not separator:
        raise ValueError("Comparison must use LATER_vs_EARLIER round names")
    combined = []
    for name in (later, earlier):
        data = load_frequency(path, name)
        wild = data.loc[data.wt_aa == data.mut_aa, ["reference_aa_position", "wt_aa", "count"]].rename(columns={"count": "wt_count"})
        data = data.loc[data.wt_aa != data.mut_aa].merge(wild, on=["reference_aa_position", "wt_aa"], how="left", validate="many_to_one")
        data["ratio"] = (data["count"] + pseudocount) / (data.wt_count.fillna(0) + pseudocount)
        combined.append(data[KEYS + ["ratio"]])
    merged = combined[0].merge(combined[1], on=KEYS, how="outer", suffixes=("_later", "_earlier"))
    # Absent mutation rows have count zero, but retain the corresponding WT denominator.
    for name, suffix in ((later, "later"), (earlier, "earlier")):
        round_data = load_frequency(path, name)
        wild = round_data.loc[round_data.wt_aa == round_data.mut_aa, ["reference_aa_position", "wt_aa", "count"]].rename(columns={"count": "wt_missing"})
        merged = merged.merge(wild, on=["reference_aa_position", "wt_aa"], how="left", validate="many_to_one")
        merged[f"ratio_{suffix}"] = merged[f"ratio_{suffix}"].fillna(pseudocount / (merged.wt_missing.fillna(0) + pseudocount))
        merged = merged.drop(columns="wt_missing")
    merged["log2_enrichment"] = np.log2(merged.ratio_later / merged.ratio_earlier)
    return merged[KEYS + ["log2_enrichment"]]


def plot_correlation(path1: Path, path2: Path, comparison: str, destination: Path, table: Path) -> None:
    frames = [enrichment_from_counts(path, comparison, 0.5) for path in (path1, path2)]
    paired = frames[0].merge(frames[1], on=KEYS, suffixes=("_rep1", "_rep2"), validate="one_to_one")
    paired = paired.replace([np.inf, -np.inf], np.nan).dropna()
    if len(paired) < 2:
        raise ValueError("At least two shared finite mutation values are required")
    paired.to_csv(table, index=False)
    x, y = paired.log2_enrichment_rep1, paired.log2_enrichment_rep2
    pearson = pearsonr(x, y).statistic; spearman = spearmanr(x, y).statistic
    fig, ax = plt.subplots(figsize=(6, 6)); ax.scatter(x, y, s=12, alpha=.5); ax.axhline(0, color='gray'); ax.axvline(0, color='gray')
    ax.set(xlabel="Replicate 1 log2 enrichment", ylabel="Replicate 2 log2 enrichment", title=f"n={len(paired)}; Pearson r={pearson:.3f}; Spearman ρ={spearman:.3f}")
    fig.tight_layout(); fig.savefig(destination, dpi=200); plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rep1-frequency', type=Path, required=True)
    parser.add_argument('--rep2-frequency', type=Path, required=True)
    parser.add_argument('--round', default='initial_round')
    parser.add_argument('--comparison', default='first_round_vs_initial_round')
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    first = load_frequency(args.rep1_frequency, args.round); second = load_frequency(args.rep2_frequency, args.round)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    plot_frequency(first, args.output_dir / 'replicate_1_frequency_heatmap.png')
    plot_frequency(second, args.output_dir / 'replicate_2_frequency_heatmap.png')
    plot_coverage(first, second, args.output_dir / 'coverage.png', args.output_dir / 'coverage.csv')
    plot_correlation(args.rep1_frequency, args.rep2_frequency, args.comparison, args.output_dir / 'replicate_correlation.png', args.output_dir / 'replicate_correlation.csv')


if __name__ == '__main__':
    main()
