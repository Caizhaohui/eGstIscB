# TAM sequence analysis

This workflow first demultiplexes pooled high-throughput sequencing reads using barcode definitions, then finds and counts TAM sequences in each sample's reads. `code/demultiplex.sh` runs Ultraplex; `code/find_TAM.py` performs counting. Input paths, barcode definitions, target sequence and output paths can be changed for each batch.

## Requirements and input

- Python 3 (standard library only for `find_TAM.py`).
- Ultraplex, available on `PATH` or specified through `ULTRAPLEX`.
- A pooled R1 FASTQ file and a barcode definition file. The example `metadata/barcodes.csv` uses generic sample names; replace it or set `BARCODES` to the appropriate file for each batch.
- The confirmed upstream target sequence and the desired TAM length.

Run the following commands from the `01_TAM_analysis` folder. The recorded PE150 analysis used R1 only. Raw and demultiplexed reads are not included in this folder.

## 1. Demultiplex reads

Set the input and a writable destination, then run the script. The barcode file defaults to `metadata/barcodes.csv`; the optional `BARCODES` variable selects another file.

```bash
export R1_FASTQ=/path/to/pooled_R1.fastq.gz
export OUTPUT_DIR=/path/to/demultiplexed
export BARCODES=metadata/barcodes.csv
export THREADS=8
bash code/demultiplex.sh
```

Ultraplex separates reads according to the barcode definitions and writes sample-specific FASTQ files into `OUTPUT_DIR`. The script retains the original `-dbr` option; check its meaning for the installed Ultraplex version. The original Slurm job requested 8 CPUs, 32 GB and one hour; adjust scheduler settings locally. Historical software versions are unknown.

## 2. Find and count TAM sequences

Run `find_TAM.py` on each sample's demultiplexed FASTQ. For example, replace the target sequence with the confirmed sequence for the experiment:

```bash
export TARGET_SEQ=YOUR_CONFIRMED_TARGET_SEQUENCE
python code/find_TAM.py \
  --input "$OUTPUT_DIR/ultraplex_demux_sample_1.fastq.gz" \
  --output /path/to/sample_1_TAM.csv \
  --target_seq "$TARGET_SEQ" \
  --length 6
```

Repeat with each sample's actual demultiplexed filename and a separate CSV output. Input can be gzip-compressed or uncompressed FASTQ. Although the original script's help also mentions FASTA, its parser reads four lines per record and requires FASTQ.

For each read, the script finds the first exact occurrence of the uppercased target and counts the following `length` bases. It skips reads containing characters outside `ATCGN` and extracted sequences containing `N`. Read sequences are not converted to uppercase. The output columns are `Subsequence`, `Count`, and `Proportion`; proportion is the count divided by all retained counts, formatted to six decimal places. The script does not calculate enrichment between samples.

The exact target sequence and command used for the historical 525 mix result tables were not recorded in the source directory; reproduction of those tables remains unverified. The counting logic is unchanged from the original `find_pam.py`; the file was renamed and its comments were translated to English.
