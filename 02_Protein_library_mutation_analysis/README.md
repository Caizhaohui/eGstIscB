# Protein library mutation analysis: minibwa workflow

This workflow uses minibwa to monitor mutation frequencies during directed protein evolution. Paired PE150 reads are demultiplexed by dual barcodes, quality-filtered with fastp, mapped to a coding reference with minibwa and counted by nucleotide, codon and amino acid. Analysis scripts preserve the original counting rules; the fastp wrapper has configurable executable/library paths, and demultiplexing requires explicit metadata and read paths. No original server paths are required.

## Code order

Files `01`–`04` prepare and count sequencing reads. File `05` visualizes minibwa results from two replicates. Files `06`–`08` form the separate local fold-change heatmap workflow. File `09` calculates fragment coverage and automatically loads its plotting helper, file `10`; run `09`, not `10` independently. The two visualization paths use different input tables and calculation rules.

## Inputs and dependencies

- Paired, synchronized R1 and R2 FASTQ files, and the matching sample metadata CSV. `metadata/example_samples.csv` provides example records with generic sample names, including the barcode pairs, PCR products, primers, fragment/sample IDs and library rounds. Replace all experiment-specific fields for a different batch and verify each PCR product occurs exactly in the chosen reference. The counter expects rounds named `initial_round`, `first_round`, and `second_round` for its three comparisons.
- `references/reference.fasta` is an example coding reference. Use the corresponding coding reference for other experiments; reference and primers determine positions and wild-type codons.
- Python with Biopython; `pigz` for compressed input with parallel decompression; fastp; [minibwa](https://github.com/lh3/minibwa); and Python packages pandas, numpy, scipy and matplotlib for plots. One original minibwa-job log recorded Python 3.9.18, Biopython 1.85 and fastp 1.0.1; minibwa reported `0.3-r396-dirty`, whose local changes are unknown. Other historical versions are unverified. A separately built minibwa may differ.

Run from this folder using a new writable work directory. Replace the example metadata and reference with files for the experiment. Raw sequencing files, intermediates and result tables are stored separately.

## 1. Set paths and prepare reference

```bash
export R1=/path/to/pooled_R1.fastq.gz R2=/path/to/pooled_R2.fastq.gz
export WORK=/path/to/new_work_directory
export METADATA="$PWD/metadata/example_samples.csv"
export REFERENCE="$PWD/references/reference.fasta"
export MINIBWA=/path/to/minibwa
mkdir -p "$WORK/reference" "$WORK/results" "$WORK/demux" "$WORK/clean"
python code/01_prepare_reference_and_fragments.py --reference "$REFERENCE" --metadata "$METADATA" \
  --out-reference "$WORK/reference/reference.fasta" --results-dir "$WORK/results"
"$MINIBWA" index -t 8 "$WORK/reference/reference.fasta"
```

Reference preparation checks exact fragment placement and writes reference, codon, fragment and primer mapping QC tables. The index is written beside the staged reference.

## 2. Demultiplex and quality-filter

```bash
python code/02_demultiplex.py --metadata "$METADATA" --r1 "$R1" --r2 "$R2" \
  --outdir "$WORK/demux" --results-dir "$WORK/results" \
  --pigz-threads 8 --output-compression none --progress-every 1000000
FASTP_BIN=fastp bash code/03_run_fastp_qc.sh "$WORK/demux" "$WORK/clean" 4 4
```

The demultiplexer requires exact 8-base prefixes on both mates, accepts either barcode orientation and removes the prefixes and corresponding quality characters. It writes paired sample files and barcode QC tables; unmatched pairs are counted. It assumes equal, synchronized input files and does not verify mate IDs. `--max-reads` defaults to zero (all reads); a positive value takes the first N pairs, not a random sample. fastp runs paired-read QC with four threads per job and four concurrent jobs, using the installed version's remaining defaults. Existing nonempty cleaned pairs are skipped: use a fresh work directory to avoid stale outputs.

## 3. Map and count mutations

```bash
python code/04_count_codon_from_minibwa_stream.py \
  --reference "$WORK/reference/reference.fasta" --metadata "$METADATA" \
  --clean-dir "$WORK/clean" --results-dir "$WORK/results" \
  --minibwa-bin "$MINIBWA" --workers 8 --map-threads 2 --pseudocount 0.5
```

The counter runs minibwa on each cleaned pair and reads its SAM stream; it does not consume the separately generated BAM files from the historical job. Bases from the two mates are joined by query name; discordant overlapping bases become `N`, and incomplete or ambiguous codons are excluded. Positions in count tables are 1-based reference positions. PCR fragments can overlap, so round-level counts combine observations across fragments. Per-position amino-acid frequency is count divided by depth. These counts should not be interpreted as full-length haplotypes.

Outputs include `base_counts_by_sample.csv`, `aa_counts_by_sample.csv`, `stream_codon_count_qc.csv`, `aa_frequency_by_library_round.csv`, `top_mutations_by_library_round.csv` and `round_enrichment_summary.csv`. The enrichment table compares first/initial, second/first and second/initial rounds using `log2(((mut_later + 0.5)/(WT_later + 0.5))/((mut_earlier + 0.5)/(WT_earlier + 0.5)))`; 0.5 is a count pseudocount. **Limitation:** the enrichment code indexes rows by amino acid and can overwrite synonymous-codon rows instead of summing them. Preserve this behavior when comparing historical outputs and evaluate its effect before interpreting enrichment.

## 4. Visualize two biological replicates

Run the analysis above independently for each biological replicate, keeping separate work directories and preserving matching reference coordinates and round definitions. Then use the two `aa_frequency_by_library_round.csv` files:

```bash
python code/05_plot_frequency_replicate_and_coverage.py \
  --rep1-frequency /path/to/replicate_1/results/aa_frequency_by_library_round.csv \
  --rep2-frequency /path/to/replicate_2/results/aa_frequency_by_library_round.csv \
  --round initial_round --comparison first_round_vs_initial_round \
  --output-dir /path/to/figures
```

This produces per-replicate mutation-frequency heatmaps, a coverage plot and CSV (number of distinct non-WT amino acids at each reference position detected at counts ≥1 and ≥5), and a replicate correlation plot and paired CSV. The heatmaps and coverage aggregate synonymous-codon counts and frequencies by amino acid. Missing position/amino-acid combinations in heatmaps are blank; coverage includes positions present in the selected input rows. Correlation uses mutations shared between the two independently processed biological replicates and reports Pearson and Spearman coefficients. It recomputes WT-normalized log2 enrichment from synonymous-codon-summed counts with a 0.5 count pseudocount, rather than using the historical enrichment table that may overwrite synonymous codons. This produces minibwa based frequency plots. The local fragment workflow instead plots log2 fold change and applies coverage thresholds per fragment before taking their union. The `--round` controls heatmap/coverage and `--comparison` specifies later and earlier rounds for correlation. Use the same reference and coordinates for both replicates.

The counter's `plot_aa_position` is experiment-specific: it hides reference aa2–8 and aa426–432, subtracts 7 from aa9 onward and excludes the terminal stop. The new visualizations use `reference_aa_position` instead; adapt any coordinate conversion explicitly for a different reference. The analysis does not infer a new background reference.

## Scope and validation

The archived count script runs minibwa on cleaned pairs and reads the SAM stream directly. The historical workflow also generated sorted BAM files. Paths are supplied as command-line arguments or environment variables.

## Local plotting workflows: log2 fold change and fragment coverage

These workflows preserve the local analysis definitions and accept different inputs from the minibwa frequency plotting script above.

### Log2 fold-change heatmap

`code/06_compare_fragment_counts.R` extracts the local per-fragment comparison: it retains counts ≥5 in each input, normalizes each retained count by its own table total, inner-joins mutations, computes log2 frequency ratio with a 1e-6 frequency pseudocount, and selects the largest absolute value for duplicate mutation labels. Run `Rscript code/06_compare_fragment_counts.R /path/to/earlier.csv /path/to/later.csv /path/to/fragment_log2fc.csv`. Prepare one table per fragment, then merge them in fragment order using `code/07_merge_fragment_foldchanges.R`. For each mutation it keeps the row with the greatest absolute log2 fold change; ties keep the first input row. The resulting `Fragment` records the input order. The original results include 2-vs-0 and 3-vs-0 merged tables; the 3-vs-2 figure input is absent from the inspected local folder, although all eight 3-vs-2 fragment tables and rendered figures are present.

Merge fragment results in the desired order (expand the input list for all fragments):

```bash
Rscript code/07_merge_fragment_foldchanges.R /path/to/merged.csv /path/to/fragment1.csv /path/to/fragment2.csv
```

`code/08_plot_log2fc_heatmap.R` is adapted from the local optimized heatmap. It reads a CSV with `Mutation` (e.g. `A42G`) and `log2_foldchange`, plus a **protein** reference FASTA, and writes PDF and PNG files:

```bash
Rscript code/08_plot_log2fc_heatmap.R /path/to/log2_foldchange.csv /path/to/protein_reference.fasta /path/to/figures comparison
```

Install R packages `dplyr`, `readr`, `stringr`, `tidyr`, `tibble`, `pheatmap`, and `gtable`, plus Cairo/Arial support for the original rendering. The retained plot design assumes mutation positions 2–424 split at 213/214, a historical reference lookup offset of +6, a stop symbol `X`, and the original colour scale. These are experiment-specific settings: review them against the input and reference before use. **Coordinate check:** the original heatmap default of +6 is retained for historical fidelity, but in the inspected example +7 matched all 5,658 mutation WT labels while +6 matched 482 of 5,657 in-range rows. Set `REFERENCE_OFFSET=7` before `Rscript` when that mapping is appropriate, and verify WT letters against your own reference. This alters WT masking compared with the historical plot. This script visualizes existing fold changes; it does not compute them from minibwa tables.

### Fragment-level amino-acid coverage

`code/09_calculate_fragment_coverage.py` and `code/10_plot_fragment_coverage.py` retain the local coverage procedure. Prepare a working directory with `referenceAA.fasta` (one protein sequence) and eight `raw/sample-0-{1..8}_AAcounts.csv` files, each containing `Mutation,Count`. Run:

```bash
uv run code/09_calculate_fragment_coverage.py --root /path/to/coverage_work
```

Outputs are written under the working directory in `results/library_coverage/` and `figures/library_coverage/`. The original calculation maps data position to FASTA position with +7, excludes FASTA positions 1–8, checks wild-type letters against the reference, and applies count thresholds 1 and 5 **to each fragment before taking their union**. It records maximum and summed fragment counts separately; summed counts are not used to decide threshold detection. These fixed reference-coordinate and eight-fragment assumptions must match the experiment. The script requires Python ≥3.13 and the dependencies listed in its inline metadata.

The local coverage and heatmap workflows use fragment `Mutation,Count` tables or merged `Mutation,log2_foldchange` tables, respectively. These differ from minibwa output tables.
