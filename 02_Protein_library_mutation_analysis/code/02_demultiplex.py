#!/usr/bin/env python3
import argparse
import csv
import gzip
import os
import subprocess
from collections import Counter


def parse_args():
    parser = argparse.ArgumentParser(description="Demultiplex paired FASTQ by dual barcodes.")
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--r1", required=True)
    parser.add_argument("--r2", required=True)
    parser.add_argument("--outdir", default="demux")
    parser.add_argument("--results-dir", default="results")
    parser.add_argument("--barcode-length", type=int, default=8)
    parser.add_argument("--progress-every", type=int, default=1000000)
    parser.add_argument("--max-reads", type=int, default=0)
    parser.add_argument("--pigz-threads", type=int, default=4)
    parser.add_argument("--output-compression", choices=["none", "gzip"], default="none")
    return parser.parse_args()


def load_metadata(path):
    rows = list(csv.DictReader(open(path, newline="")))
    lookup = {}
    for row in rows:
        sample = row["SampleNumber"]
        b1 = row["Barcode1"].upper().encode()
        b2 = row["Barcode2"].upper().encode()
        lookup[(b2, b1)] = (sample, "R1_Barcode2_R2_Barcode1")
        lookup[(b1, b2)] = (sample, "R1_Barcode1_R2_Barcode2")
    return rows, lookup


def open_fastq_stream(path, pigz_threads):
    if path.endswith(".gz") and pigz_threads > 0:
        proc = subprocess.Popen(
            ["pigz", "-dc", "-p", str(pigz_threads), path],
            stdout=subprocess.PIPE,
            bufsize=1024 * 1024,
        )
        return proc.stdout, proc
    if path.endswith(".gz"):
        return gzip.open(path, "rb"), None
    return open(path, "rb", buffering=1024 * 1024), None


def fastq_iter(path, pigz_threads):
    handle, proc = open_fastq_stream(path, pigz_threads)
    try:
        while True:
            title = handle.readline()
            if not title:
                break
            seq = handle.readline()
            plus = handle.readline()
            qual = handle.readline()
            if not qual:
                break
            yield title, seq, plus, qual
    finally:
        handle.close()
        if proc is not None:
            proc.wait()


def open_output(path, compression):
    if compression == "gzip":
        return gzip.open(path + ".gz", "wb", compresslevel=1)
    return open(path, "wb", buffering=1024 * 1024)


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    os.makedirs(args.results_dir, exist_ok=True)

    rows, lookup = load_metadata(args.metadata)
    sample_handles = {}
    sample_counts = Counter()
    orient_counts = Counter()
    unmatched_counts = Counter()
    total = 0
    matched = 0

    try:
        for r1_rec, r2_rec in zip(
            fastq_iter(args.r1, args.pigz_threads),
            fastq_iter(args.r2, args.pigz_threads),
        ):
            if args.max_reads and total >= args.max_reads:
                break
            total += 1
            r1_title, r1_seq, r1_plus, r1_qual = r1_rec
            r2_title, r2_seq, r2_plus, r2_qual = r2_rec
            r1_barcode = r1_seq[:args.barcode_length].upper()
            r2_barcode = r2_seq[:args.barcode_length].upper()
            hit = lookup.get((r1_barcode, r2_barcode))
            if hit is None:
                unmatched_counts[(r1_barcode, r2_barcode)] += 1
                continue

            sample, orientation = hit
            matched += 1
            orient_counts[orientation] += 1
            sample_counts[sample] += 1

            if sample not in sample_handles:
                r1_path = os.path.join(args.outdir, f"{sample}_R1.fastq")
                r2_path = os.path.join(args.outdir, f"{sample}_R2.fastq")
                sample_handles[sample] = (
                    open_output(r1_path, args.output_compression),
                    open_output(r2_path, args.output_compression),
                )
            out1, out2 = sample_handles[sample]

            if orientation == "R1_Barcode2_R2_Barcode1":
                s1_title, s1_seq, s1_plus, s1_qual = r1_title, r1_seq, r1_plus, r1_qual
                s2_title, s2_seq, s2_plus, s2_qual = r2_title, r2_seq, r2_plus, r2_qual
            else:
                s1_title, s1_seq, s1_plus, s1_qual = r2_title, r2_seq, r2_plus, r2_qual
                s2_title, s2_seq, s2_plus, s2_qual = r1_title, r1_seq, r1_plus, r1_qual

            out1.write(s1_title)
            out1.write(s1_seq[args.barcode_length:])
            out1.write(s1_plus)
            out1.write(s1_qual[args.barcode_length:])

            out2.write(s2_title)
            out2.write(s2_seq[args.barcode_length:])
            out2.write(s2_plus)
            out2.write(s2_qual[args.barcode_length:])

            if args.progress_every and total % args.progress_every == 0:
                print(f"processed_pairs={total} matched={matched} unmatched={total-matched}", flush=True)
    finally:
        for handle1, handle2 in sample_handles.values():
            handle1.close()
            handle2.close()

    with open(os.path.join(args.results_dir, "barcode_demux_qc.csv"), "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["metric", "value"])
        writer.writeheader()
        writer.writerows([
            {"metric": "total_pairs", "value": total},
            {"metric": "matched_pairs", "value": matched},
            {"metric": "unmatched_pairs", "value": total - matched},
            {"metric": "match_rate", "value": matched / total if total else 0},
        ])

    with open(os.path.join(args.results_dir, "barcode_orientation_qc.csv"), "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["orientation", "read_pairs"])
        writer.writeheader()
        for key, value in sorted(orient_counts.items()):
            writer.writerow({"orientation": key, "read_pairs": value})

    with open(os.path.join(args.results_dir, "sample_read_counts.csv"), "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["SampleNumber", "library_round", "fragment_id", "read_pairs"])
        writer.writeheader()
        for row in rows:
            writer.writerow({
                "SampleNumber": row["SampleNumber"],
                "library_round": row["library_round"],
                "fragment_id": row["SampleNumber"].split("-")[-1],
                "read_pairs": sample_counts[row["SampleNumber"]],
            })

    with open(os.path.join(args.results_dir, "unmatched_barcode_counts.csv"), "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["R1_barcode", "R2_barcode", "count"])
        writer.writeheader()
        for (b1, b2), count in unmatched_counts.most_common(1000):
            writer.writerow({"R1_barcode": b1.decode(), "R2_barcode": b2.decode(), "count": count})


if __name__ == "__main__":
    main()
