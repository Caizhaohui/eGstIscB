#!/usr/bin/env python3
import argparse
import csv
import gzip
import math
import os
import subprocess
import sys
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

from Bio import SeqIO


CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}
VALID_BASES = set("ACGT")
CIGAR_OPS = "MIDNSHP=XB"


def normalize(seq):
    return "".join(seq.split()).upper()


def write_csv(path, fieldnames, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def load_reference(path):
    record = next(SeqIO.parse(path, "fasta"))
    seq = normalize(str(record.seq))
    codons = [seq[i:i + 3] for i in range(0, len(seq), 3)]
    aas = [CODON_TABLE.get(c, "X") for c in codons]
    return record.id, seq, codons, aas


def load_metadata(path):
    return list(csv.DictReader(open(path, newline="")))


def add_fragment_coordinates(rows, ref_seq):
    out = []
    for row in rows:
        item = dict(row)
        pcr_product = normalize(row["PCR_product"])
        start0 = ref_seq.find(pcr_product)
        if start0 < 0:
            raise ValueError(f"{row['SampleNumber']} PCR_product does not map exactly to reference")
        item["fragment_start0"] = start0
        item["fragment_end0"] = start0 + len(pcr_product)
        out.append(item)
    return out


def parse_cigar(cigar):
    num = []
    for char in cigar:
        if char.isdigit():
            num.append(char)
        elif char in CIGAR_OPS:
            if not num:
                raise ValueError(f"bad CIGAR: {cigar}")
            yield int("".join(num)), char
            num = []
    if num:
        raise ValueError(f"bad CIGAR tail: {cigar}")


def sam_record_to_ref_bases(fields, fragment_start0, fragment_end0):
    flag = int(fields[1])
    if flag & 0x4 or flag & 0x100 or flag & 0x800:
        return None
    pos0 = int(fields[3]) - 1
    cigar = fields[5]
    seq = fields[9].upper()
    if cigar == "*" or seq == "*":
        return None
    ref_pos = pos0
    query_pos = 0
    bases = {}
    for length, op in parse_cigar(cigar):
        if op in ("M", "=", "X"):
            for i in range(length):
                base = seq[query_pos + i]
                if base in VALID_BASES and fragment_start0 <= ref_pos + i < fragment_end0:
                    bases[ref_pos + i] = base
            ref_pos += length
            query_pos += length
        elif op == "I":
            query_pos += length
        elif op in ("D", "N"):
            ref_pos += length
        elif op == "S":
            query_pos += length
        elif op in ("H", "P"):
            continue
    return bases


def add_pair_counts(read_bases, codon_count, codon_counts, base_counts, metrics):
    if not read_bases:
        return
    metrics["read_pairs_seen"] += 1
    metrics[f"read_groups_with_{min(len(read_bases), 3)}_mapped_records"] += 1
    merged = {}
    discordant_bases = 0
    for bases in read_bases:
        for rpos, base in bases.items():
            old = merged.get(rpos)
            if old is None:
                merged[rpos] = base
            elif old != base:
                discordant_bases += 1
                merged[rpos] = "N"
    if discordant_bases:
        metrics["discordant_mate_bases"] += discordant_bases
    valid_positions = [pos for pos, base in merged.items() if base in VALID_BASES]
    if not valid_positions:
        return
    for rpos in valid_positions:
        base_counts[(rpos + 1, merged[rpos])] += 1
    first_codon = max(0, (min(valid_positions) + 2) // 3)
    last_codon = min(codon_count - 1, max(valid_positions) // 3)
    for codon_idx in range(first_codon, last_codon + 1):
        start = codon_idx * 3
        try:
            codon = merged[start] + merged[start + 1] + merged[start + 2]
        except KeyError:
            continue
        if all(base in VALID_BASES for base in codon):
            codon_counts[(codon_idx + 1, codon)] += 1
            metrics["codon_observations_counted"] += 1
        else:
            metrics["discordant_mate_codon_observations"] += 1


def run_minibwa_sample(task):
    row, reference, clean_dir, minibwa_bin, threads, wt_codons, wt_aas, ref_seq = task
    sample = row["SampleNumber"]
    r1 = os.path.join(clean_dir, f"{sample}_R1.clean.fastq.gz")
    r2 = os.path.join(clean_dir, f"{sample}_R2.clean.fastq.gz")
    if not os.path.exists(r1) or not os.path.exists(r2):
        raise FileNotFoundError(f"missing clean FASTQ for {sample}")
    cmd = [minibwa_bin, "map", "-t", str(threads), reference, r1, r2]
    log_dir = os.path.join("analysis_minibwa_bam", "logs", "minibwa-stream")
    os.makedirs(log_dir, exist_ok=True)
    stderr_path = os.path.join(log_dir, f"{sample}.minibwa.stderr.log")
    stderr_handle = open(stderr_path, "w")
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=stderr_handle, text=True, bufsize=1024 * 1024)
    codon_counts = Counter()
    base_counts = Counter()
    metrics = Counter()
    current_qname = None
    current_bases = []

    for line in proc.stdout:
        if not line or line[0] == "@":
            continue
        fields = line.rstrip("\n").split("\t")
        metrics["sam_records"] += 1
        qname = fields[0]
        if current_qname is not None and qname != current_qname:
            add_pair_counts(current_bases, len(wt_codons), codon_counts, base_counts, metrics)
            current_bases = []
        current_qname = qname
        ref_bases = sam_record_to_ref_bases(fields, row["fragment_start0"], row["fragment_end0"])
        if ref_bases is None:
            metrics["unmapped_or_filtered_records"] += 1
            continue
        metrics["mapped_primary_records"] += 1
        if not ref_bases:
            metrics["mapped_records_without_fragment_bases"] += 1
        current_bases.append(ref_bases)

    if current_qname is not None:
        add_pair_counts(current_bases, len(wt_codons), codon_counts, base_counts, metrics)

    return_code = proc.wait()
    stderr_handle.close()
    if return_code != 0:
        tail = ""
        try:
            with open(stderr_path) as handle:
                tail = "".join(handle.readlines()[-20:])
        except OSError:
            pass
        raise RuntimeError(f"minibwa failed for {sample} with code {return_code}\n{tail}")

    sample_total_depth = Counter()
    for (aa_pos, codon), count in codon_counts.items():
        sample_total_depth[aa_pos] += count
    aa_rows = []
    for (aa_pos, codon), count in sorted(codon_counts.items()):
        wt_codon = wt_codons[aa_pos - 1]
        wt_aa = wt_aas[aa_pos - 1]
        mut_aa = CODON_TABLE.get(codon, "X")
        depth = sample_total_depth[aa_pos]
        aa_rows.append({
            "SampleNumber": sample,
            "library_round": row["library_round"],
            "fragment_id": sample.split("-")[-1],
            "aa_position": aa_pos,
            "wt_codon": wt_codon,
            "mut_codon": codon,
            "wt_aa": wt_aa,
            "mut_aa": mut_aa,
            "count": count,
            "depth": depth,
            "frequency": count / depth if depth else 0.0,
            "is_mutation": codon != wt_codon,
            "aa_mutation": f"{wt_aa}{aa_pos}{mut_aa}",
        })

    nt_depth = Counter()
    for (nt_pos, base), count in base_counts.items():
        nt_depth[nt_pos] += count
    base_rows = []
    for (nt_pos, base), count in sorted(base_counts.items()):
        ref_base = ref_seq[nt_pos - 1]
        depth = nt_depth[nt_pos]
        base_rows.append({
            "SampleNumber": sample,
            "library_round": row["library_round"],
            "fragment_id": sample.split("-")[-1],
            "nt_position": nt_pos,
            "ref_base": ref_base,
            "base": base,
            "count": count,
            "depth": depth,
            "frequency": count / depth if depth else 0.0,
        })

    metric_row = {
        "SampleNumber": sample,
        "library_round": row["library_round"],
        "fragment_id": sample.split("-")[-1],
        "minibwa_return_code": return_code,
        "minibwa_stderr_log": stderr_path,
    }
    metric_row.update(metrics)
    return sample, aa_rows, base_rows, metric_row


def plot_position(aa_pos, wt_aa):
    if aa_pos == 1:
        return 1
    if 2 <= aa_pos <= 8 or 426 <= aa_pos <= 432 or wt_aa == "*":
        return ""
    return aa_pos - 7


def summarize_rounds(sample_rows):
    grouped = Counter()
    depths = Counter()
    mut_codons = defaultdict(set)
    for row in sample_rows:
        aa_pos = int(row["aa_position"])
        key = (row["library_round"], aa_pos, row["wt_codon"], row["wt_aa"], row["mut_aa"])
        grouped[key] += int(row["count"])
        mut_codons[key].add(row["mut_codon"])
        depths[(row["library_round"], aa_pos)] += int(row["count"])
    rows = []
    for (round_name, aa_pos, wt_codon, wt_aa, mut_aa), count in sorted(grouped.items()):
        depth = depths[(round_name, aa_pos)]
        rows.append({
            "library_round": round_name,
            "aa_position": aa_pos,
            "reference_aa_position": aa_pos,
            "plot_aa_position": plot_position(aa_pos, wt_aa),
            "wt_codon": wt_codon,
            "mut_codon": ";".join(sorted(mut_codons[(round_name, aa_pos, wt_codon, wt_aa, mut_aa)])),
            "wt_aa": wt_aa,
            "mut_aa": mut_aa,
            "count": count,
            "depth": depth,
            "frequency": count / depth if depth else 0.0,
            "is_mutation": mut_aa != wt_aa,
            "aa_mutation": f"{wt_aa}{aa_pos}{mut_aa}",
        })
    return rows


def summarize_top(round_rows):
    rows = [row for row in round_rows if row["is_mutation"] is True or str(row["is_mutation"]).lower() == "true"]
    rows.sort(key=lambda r: (r["library_round"], -float(r["frequency"]), int(r["aa_position"]), r["mut_aa"]))
    return rows


def enrichment_rows(round_rows, pseudocount):
    by_key = {}
    depth_by_round_pos = {}
    wt_by_round_pos = {}
    for row in round_rows:
        round_name = row["library_round"]
        aa_pos = int(row["aa_position"])
        key = (round_name, aa_pos, row["mut_aa"])
        by_key[key] = row
        depth_by_round_pos[(round_name, aa_pos)] = int(row["depth"])
        if row["mut_aa"] == row["wt_aa"]:
            wt_by_round_pos[(round_name, aa_pos)] = row
    comparisons = [("first_round", "initial_round"), ("second_round", "first_round"), ("second_round", "initial_round")]
    out = []
    for num_round, den_round in comparisons:
        keys = {(int(r["aa_position"]), r["mut_aa"]) for r in round_rows if r["library_round"] in {num_round, den_round}}
        for aa_pos, mut_aa in keys:
            num = by_key.get((num_round, aa_pos, mut_aa))
            den = by_key.get((den_round, aa_pos, mut_aa))
            template = num or den
            if template["mut_aa"] == template["wt_aa"]:
                continue
            num_wt = wt_by_round_pos.get((num_round, aa_pos))
            den_wt = wt_by_round_pos.get((den_round, aa_pos))
            num_count = int(num["count"]) if num else 0
            den_count = int(den["count"]) if den else 0
            num_wt_count = int(num_wt["count"]) if num_wt else 0
            den_wt_count = int(den_wt["count"]) if den_wt else 0
            num_depth = depth_by_round_pos.get((num_round, aa_pos), 0)
            den_depth = depth_by_round_pos.get((den_round, aa_pos), 0)
            num_freq = num_count / num_depth if num_depth else 0.0
            den_freq = den_count / den_depth if den_depth else 0.0
            num_wt_freq = num_wt_count / num_depth if num_depth else 0.0
            den_wt_freq = den_wt_count / den_depth if den_depth else 0.0
            num_mut_vs_wt = (num_count + pseudocount) / (num_wt_count + pseudocount)
            den_mut_vs_wt = (den_count + pseudocount) / (den_wt_count + pseudocount)
            out.append({
                "comparison": f"{num_round}_vs_{den_round}",
                "numerator_round": num_round,
                "denominator_round": den_round,
                "aa_position": aa_pos,
                "reference_aa_position": aa_pos,
                "plot_aa_position": template["plot_aa_position"],
                "wt_aa": template["wt_aa"],
                "mut_aa": template["mut_aa"],
                "mut_codon": template["mut_codon"],
                "aa_mutation": template["aa_mutation"],
                "numerator_count": num_count,
                "numerator_depth": num_depth,
                "numerator_wt_count": num_wt_count,
                "denominator_count": den_count,
                "denominator_depth": den_depth,
                "denominator_wt_count": den_wt_count,
                "numerator_frequency": num_freq,
                "denominator_frequency": den_freq,
                "numerator_wt_frequency": num_wt_freq,
                "denominator_wt_frequency": den_wt_freq,
                "numerator_mut_vs_wt_ratio": num_mut_vs_wt,
                "denominator_mut_vs_wt_ratio": den_mut_vs_wt,
                "log2_enrichment": math.log2(num_mut_vs_wt / den_mut_vs_wt) if den_mut_vs_wt else "",
                "enrichment_method": "WT_normalized_log2_mut_vs_wt_odds_ratio",
                "is_mutation": template["is_mutation"],
            })
    return out


def main():
    parser = argparse.ArgumentParser(description="Run minibwa and count read-backed codons directly from SAM stream.")
    parser.add_argument("--reference", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--clean-dir", default="analysis_minibwa_bam/clean")
    parser.add_argument("--results-dir", default="analysis_minibwa_bam/results")
    parser.add_argument("--minibwa-bin", default="tools/minibwa/bin/minibwa")
    parser.add_argument("--map-threads", type=int, default=2)
    parser.add_argument("--workers", type=int, default=8)
    parser.add_argument("--pseudocount", type=float, default=0.5)
    args = parser.parse_args()

    os.makedirs(args.results_dir, exist_ok=True)
    _, ref_seq, wt_codons, wt_aas = load_reference(args.reference)
    rows = add_fragment_coordinates(load_metadata(args.metadata), ref_seq)
    tasks = [(row, args.reference, args.clean_dir, args.minibwa_bin, args.map_threads, wt_codons, wt_aas, ref_seq) for row in rows]

    sample_rows = []
    base_rows = []
    qc_rows = []
    with ProcessPoolExecutor(max_workers=max(1, args.workers)) as executor:
        future_to_sample = {executor.submit(run_minibwa_sample, task): task[0]["SampleNumber"] for task in tasks}
        for future in as_completed(future_to_sample):
            sample = future_to_sample[future]
            sample_name, aa_rows, nt_rows, metric_row = future.result()
            sample_rows.extend(aa_rows)
            base_rows.extend(nt_rows)
            qc_rows.append(metric_row)
            print(f"stream_counted {sample_name}: aa_rows={len(aa_rows)} base_rows={len(nt_rows)}", flush=True)

    sample_rows.sort(key=lambda r: (r["SampleNumber"], int(r["aa_position"]), r["mut_codon"]))
    base_rows.sort(key=lambda r: (r["SampleNumber"], int(r["nt_position"]), r["base"]))
    qc_rows.sort(key=lambda r: r["SampleNumber"])

    aa_fields = ["SampleNumber", "library_round", "fragment_id", "aa_position", "wt_codon", "mut_codon", "wt_aa", "mut_aa", "count", "depth", "frequency", "is_mutation", "aa_mutation"]
    base_fields = ["SampleNumber", "library_round", "fragment_id", "nt_position", "ref_base", "base", "count", "depth", "frequency"]
    round_fields = ["library_round", "aa_position", "reference_aa_position", "plot_aa_position", "wt_codon", "mut_codon", "wt_aa", "mut_aa", "count", "depth", "frequency", "is_mutation", "aa_mutation"]
    write_csv(os.path.join(args.results_dir, "aa_counts_by_sample.csv"), aa_fields, sample_rows)
    write_csv(os.path.join(args.results_dir, "base_counts_by_sample.csv"), base_fields, base_rows)
    write_csv(os.path.join(args.results_dir, "stream_codon_count_qc.csv"), sorted({key for row in qc_rows for key in row}), qc_rows)
    round_rows = summarize_rounds(sample_rows)
    write_csv(os.path.join(args.results_dir, "aa_frequency_by_library_round.csv"), round_fields, round_rows)
    write_csv(os.path.join(args.results_dir, "top_mutations_by_library_round.csv"), round_fields, summarize_top(round_rows))
    enrich = enrichment_rows(round_rows, args.pseudocount)
    write_csv(os.path.join(args.results_dir, "round_enrichment_summary.csv"), list(enrich[0]) if enrich else [], enrich)


if __name__ == "__main__":
    main()
