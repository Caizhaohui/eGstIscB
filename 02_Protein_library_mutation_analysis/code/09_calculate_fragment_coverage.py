#!/usr/bin/env -S uv run --script
# /// script
# requires-python = ">=3.13"
# dependencies = ["numpy==2.2.6", "polars==2.0.0", "pyarrow==23.0.1", "pydantic>=2.12,<3", "typer>=0.20,<1", "matplotlib==3.10.7"]
# ///
# ─── How to run ───
# 1. Install uv: https://docs.astral.sh/uv/getting-started/installation/
# 2. Run: python -m uv run code/09_calculate_fragment_coverage.py
#    Optional: --root PATH (defaults to this script's project directory)
# ──────────────────
"""Count group-0 amino-acid coverage without pooling threshold evidence."""

from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from time import perf_counter
from typing import Annotated, ClassVar, Final, override

import numpy as np
import polars as pl
import typer
from pydantic import BaseModel, ConfigDict, Field, TypeAdapter

AA: Final = "ACDEFGHIKLMNPQRSTVWY"
SYMBOLS: Final = AA + "X"
OFFSET: Final = 7
EXCLUDED_NTERM: Final = 8
THRESHOLDS: Final = (1, 5)
DEFAULT_ROOT: Final = Path(__file__).resolve().parents[1]


@dataclass(frozen=True, slots=True)
class InputError(Exception):
    """Describe a failed input contract before writing results."""

    detail: str

    @override
    def __str__(self) -> str:
        return self.detail


class MutationRow(BaseModel):
    """Parse the two columns of an original AAcounts record."""

    model_config: ClassVar[ConfigDict] = ConfigDict(frozen=True, strict=True, extra="forbid")
    Mutation: Annotated[str, Field(pattern=r"^[ACDEFGHIKLMNPQRSTVWY]-?\d+[ACDEFGHIKLMNPQRSTVWYX]$")]
    Count: Annotated[int, Field(ge=0)]


def load_inputs(root: Path) -> tuple[str, pl.DataFrame]:
    """Read exactly eight group-0 files and verify the +7 reference mapping."""
    lines = root.joinpath("referenceAA.fasta").read_text(encoding="utf-8-sig").splitlines()
    if sum(line.startswith(">") for line in lines) != 1:
        raise InputError("Reference FASTA must contain exactly one sequence.")
    seq = "".join(line.strip() for line in lines if not line.startswith(">")).upper()
    if not seq or set(seq) - set(AA):
        raise InputError("Reference must contain only the 20 standard amino acids.")
    frames: list[pl.DataFrame] = []
    adapter = TypeAdapter(tuple[MutationRow, ...])
    for fragment in range(1, 9):
        path = root / "raw" / f"sample-0-{fragment}_AAcounts.csv"
        records = adapter.validate_python(pl.read_csv(path, schema={"Mutation": pl.String, "Count": pl.Int64}).to_dicts())
        if not records or len({r.Mutation for r in records}) != len(records):
            raise InputError(f"Empty file or duplicate mutation labels: {path.name}")
        positions = [int(r.Mutation[1:-1]) + OFFSET for r in records]
        if any(not 1 <= p <= len(seq) or seq[p - 1] != r.Mutation[0] for p, r in zip(positions, records, strict=True)):
            raise InputError(f"Mutation/reference mismatch using +7 offset: {path.name}")
        frames.append(pl.DataFrame({
            "Source_File": [path.name] * len(records), "Fragment": [fragment] * len(records),
            "Mutation": [r.Mutation for r in records], "Count": [r.Count for r in records],
            "FASTA_Position": positions, "Data_Position": [p - OFFSET for p in positions],
            "WT_AA": [r.Mutation[0] for r in records], "Observed_AA": [r.Mutation[-1] for r in records],
        }))
    return seq, pl.concat(frames).with_columns(
        (pl.col("Count") >= 1).alias("Pass_Count1"), (pl.col("Count") >= 5).alias("Pass_Count5"),
    )


def count_grid(raw: pl.DataFrame, seq: str) -> pl.DataFrame:
    """Keep max and summed counts distinct on the complete 432 x 21 grid."""
    grouped = raw.group_by("FASTA_Position", "Observed_AA").agg(
        pl.col("Count").max().alias("Max_Fragment_Count"),
        pl.col("Count").sum().alias("Summed_Count"),
        pl.len().alias("Fragment_Records"),
        pl.col("Source_File").sort().str.join(";").alias("Source_Files"),
    )
    return pl.DataFrame({
        "FASTA_Position": np.repeat(np.arange(1, len(seq) + 1), len(SYMBOLS)),
        "Observed_AA": list(SYMBOLS) * len(seq), "WT_AA": np.repeat(list(seq), len(SYMBOLS)),
    }).join(grouped, on=["FASTA_Position", "Observed_AA"], how="left").with_columns(
        pl.col("Max_Fragment_Count", "Summed_Count", "Fragment_Records").fill_null(0),
        pl.col("Source_Files").fill_null(""),
        (pl.col("FASTA_Position") - OFFSET).alias("Data_Position"),
        (pl.col("Observed_AA") == pl.col("WT_AA")).alias("Is_Wildtype"),
    ).sort("FASTA_Position", "Observed_AA")


def position_table(grid: pl.DataFrame, seq: str) -> pl.DataFrame:
    """Report both thresholds in one row per reference position."""
    rows: list[dict[str, str | int | float | bool]] = []
    for pos, wt in enumerate(seq, 1):
        part = grid.filter(pl.col("FASTA_Position") == pos)
        row: dict[str, str | int | float | bool] = {
            "FASTA_Position": pos, "Data_Position": pos - OFFSET, "WT_AA": wt,
            "Has_Input_Records": bool(part["Fragment_Records"].sum()),
            "In_Data_2_424": 2 <= pos - OFFSET <= 424,
        }
        for threshold in THRESHOLDS:
            detected = set(part.filter(pl.col("Max_Fragment_Count") >= threshold)["Observed_AA"].to_list())
            present = "".join(a for a in AA if a in detected)
            missing = "".join(a for a in AA if a not in detected)
            prefix = f"Count{threshold}_"
            row.update({
                prefix + "Standard_AAs": present, prefix + "Standard_AA_Count": len(present),
                prefix + "Missing_Standard_AAs": missing, prefix + "Standard_Coverage": len(present) / len(AA),
                prefix + "WT_Detected": wt in detected,
                prefix + "NonWT_AAs": present.replace(wt, ""),
                prefix + "NonWT_AA_Count": len(present.replace(wt, "")),
                prefix + "NonWT_Coverage": len(present.replace(wt, "")) / (len(AA) - 1),
                prefix + "Stop_X_Detected": "X" in detected,
            })
        rows.append(row)
    return pl.DataFrame(rows)


def summary_tables(positions: pl.DataFrame) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Use all positions in each stated scope as the coverage denominator."""
    summaries: list[dict[str, str | int | float]] = []
    distributions: list[dict[str, str | int]] = []
    scopes = (("Reference_excluding_first8", positions), ("Data_2_424", positions.filter(pl.col("In_Data_2_424"))))
    for scope, data in scopes:
        for threshold in THRESHOLDS:
            prefix = f"Count{threshold}_"
            standard = TypeAdapter(list[int]).validate_python(data[prefix + "Standard_AA_Count"].to_list())
            mutant = TypeAdapter(list[int]).validate_python(data[prefix + "NonWT_AA_Count"].to_list())
            summaries.append({
                "Scope": scope, "Count_Threshold": threshold, "Positions": data.height,
                "Positions_With_Input": data.filter(pl.col("Has_Input_Records")).height,
                "Detected_Standard_Combinations": sum(standard),
                "Possible_Standard_Combinations": data.height * len(AA),
                "Standard_Coverage": sum(standard) / (data.height * len(AA)),
                "Detected_NonWT_Combinations": sum(mutant),
                "Possible_NonWT_Combinations": data.height * (len(AA) - 1),
                "NonWT_Coverage": sum(mutant) / (data.height * (len(AA) - 1)),
                "Positions_With_All20": standard.count(len(AA)),
                "Positions_With_All19_NonWT": mutant.count(len(AA) - 1),
                "Positions_With_WT": data.filter(pl.col(prefix + "WT_Detected")).height,
                "Positions_With_Stop_X": data.filter(pl.col(prefix + "Stop_X_Detected")).height,
                "Mean_Standard_AA_Count": sum(standard) / data.height,
            })
            for kind, values, maximum in (("Standard", standard, len(AA)), ("NonWT", mutant, len(AA) - 1)):
                distributions.extend({"Scope": scope, "Count_Threshold": threshold, "Kind": kind,
                                      "AA_Count": n, "Position_Count": values.count(n)}
                                     for n in range(maximum + 1))
    return pl.DataFrame(summaries), pl.DataFrame(distributions)


def write_report(root: Path, summary: pl.DataFrame) -> None:
    """Document denominators, provenance, and important limits in Chinese."""
    sections = [
        "# 0 组蛋白突变文库覆盖度\n",
        "输入：raw/sample-0-{1..8}_AAcounts.csv；参考：referenceAA.fasta。",
        "按实验设计排除未进行突变的 N 端 8 个氨基酸 MPKKKRKV（FASTA 1–8，数据编号 −6 至 1），不纳入任何覆盖度统计、分母或图表。保留原始坐标，不重新编号。",
        "坐标：FASTA 从 1 开始，FASTA_Position = Data_Position + 7；所有输入野生型字母均已核对一致。",
        "阈值：任一原始片段 Count ≥ 阈值即检出；重叠片段取种类并集。Summed_Count 仅提供计数明细，不能用于判定 ≥5。",
        "20 种标准氨基酸包括野生型；非野生型分母为每个位点 19 种。X 按项目既有约定视为终止子，单独报告。",
        "排除前 8 位后，主要范围为 FASTA 9–432（数据编号 2–425），共 424 位，分母为 424×20（非野生型为 424×19）；另汇总数据 2–424 范围，即 FASTA 9–431，共 423 位。",
        "FASTA 第 432 位在 0 组无输入记录，保留在 424 位范围中并标记无数据；423 位范围不含该位置，分别列示以明确分母。",
        "Has_Input_Records=false 表示该位置没有输入记录，不能推断生物学上不存在。缺失种类仅表示未达到当前阈值。",
        "本统计是测序检出覆盖度，不能从这些汇总表推断独立克隆数、分子数或多位点突变组合覆盖度。\n",
    ]
    for record in summary.iter_rows(named=True):
        sections.append(f"- {record['Scope']}，Count ≥ {record['Count_Threshold']}：标准组合 {record['Detected_Standard_Combinations']}/{record['Possible_Standard_Combinations']}（{record['Standard_Coverage']:.2%}）；非野生型组合 {record['Detected_NonWT_Combinations']}/{record['Possible_NonWT_Combinations']}（{record['NonWT_Coverage']:.2%}）；20 种齐全位点 {record['Positions_With_All20']}/{record['Positions']}。")
    sections.extend([
        "\n## 文件与复现",
        "- per_position_coverage.csv：424 行（FASTA 9–432），两种阈值的逐位点结果；氨基酸列表按 ACDEFGHIKLMNPQRSTVWY 排序。",
        "- amino_acid_counts.csv：完整 424×21 网格，包含原始片段最大计数、计数和及来源。",
        "- source_mutation_counts.csv：排除 FASTA 1–8 后的逐来源记录及坐标、阈值标记；原始输入文件保持不变。",
        "- coverage_summary.csv、coverage_distribution.csv：覆盖率及种类数分布。",
        "- fragment_audit.csv：各片段范围、记录数及阈值保留记录数；这不是独立 reads 总数。",
        "- validation.txt：真实输出的检查结果；figures/library_coverage/：PNG、PDF 图表。",
        "\n在项目目录运行：`python -m uv run code/09_calculate_fragment_coverage.py`。脚本只写入本次新增输出目录，重复运行会更新这些结果。",
        "现有热图脚本使用 +6，本次统计使用经原始记录验证的 +7；现有热图及其脚本未修改。",
    ])
    _ = root.joinpath("results/library_coverage/README.md").write_text("\n\n".join(sections) + "\n", encoding="utf-8")


def main(root: Path = DEFAULT_ROOT) -> None:
    """Generate the group-0 coverage tables and publication-exportable figures."""
    start = perf_counter()
    sys.path.insert(0, str(DEFAULT_ROOT))
    from importlib import import_module
    draw_figures = import_module("10_plot_fragment_coverage").draw_figures

    seq, original = load_inputs(root)
    raw = original.filter(pl.col("FASTA_Position") > EXCLUDED_NTERM)
    grid = count_grid(raw, seq).filter(pl.col("FASTA_Position") > EXCLUDED_NTERM)
    positions = position_table(grid, seq).filter(pl.col("FASTA_Position") > EXCLUDED_NTERM)
    summary, distribution = summary_tables(positions)
    checks = [
        (positions.height == len(seq) - EXCLUDED_NTERM, "One output row per eligible FASTA position (9 to end)"),
        (grid.height == (len(seq) - EXCLUDED_NTERM) * len(SYMBOLS), "Complete eligible-position x 21-symbol grid"),
        (positions.filter(pl.col("FASTA_Position") <= EXCLUDED_NTERM).is_empty(), "First 8 FASTA positions excluded"),
        (grid.select("FASTA_Position", "Observed_AA").is_duplicated().sum() == 0, "No duplicate position/AA cells"),
        (positions.filter(pl.col("Count1_Standard_AA_Count") > 20).is_empty(), "At most 20 standard AA types"),
        (positions.filter(pl.col("Count1_NonWT_AA_Count") > 19).is_empty(), "At most 19 non-WT AA types"),
        (positions.filter(pl.col("Count5_Standard_AA_Count") > pl.col("Count1_Standard_AA_Count")).is_empty(), "Count5 counts do not exceed Count1"),
    ]
    if any(not passed for passed, _ in checks):
        raise InputError("Output invariant failed: " + "; ".join(label for passed, label in checks if not passed))
    output = root / "results/library_coverage"
    output.mkdir(parents=True, exist_ok=True)
    for name, table in (("per_position_coverage", positions), ("amino_acid_counts", grid),
                        ("source_mutation_counts", raw), ("coverage_summary", summary),
                        ("coverage_distribution", distribution)):
        table.write_csv(output / f"{name}.csv", include_bom=True)
    raw.group_by("Source_File", "Fragment").agg(
        pl.len().alias("Records"), pl.col("Data_Position").min().alias("Min_Data_Position"),
        pl.col("Data_Position").max().alias("Max_Data_Position"),
        pl.col("FASTA_Position").n_unique().alias("Positions"),
        pl.col("Pass_Count1", "Pass_Count5").sum(),
    ).sort("Fragment").write_csv(output / "fragment_audit.csv", include_bom=True)
    write_report(root, summary)
    draw_figures(root, positions, grid)
    _ = (output / "validation.txt").write_text(f"Input checks passed: {original.height} records in 8 group-0 files; valid counts and labels; no duplicate labels; all WT letters match FASTA with +7 offset.\nExcluded {original.height - raw.height} records at FASTA positions 1-8 before statistics.\n" + "\n".join(f"PASS: {label}" for _, label in checks) + "\n", encoding="utf-8")
    typer.echo(summary.write_csv())
    typer.echo(f"Processed {raw.height:,} eligible records; {positions.height} positions (FASTA 9-{len(seq)}); {perf_counter() - start:.2f} seconds.")


if __name__ == "__main__":
    typer.run(main)
