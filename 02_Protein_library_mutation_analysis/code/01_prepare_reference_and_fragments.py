#!/usr/bin/env python3
import argparse
import csv
import os

from Bio import SeqIO
from Bio.Seq import Seq


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


def normalize(seq):
    return "".join(seq.split()).upper()


def write_csv(path, fieldnames, rows):
    with open(path, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description="Prepare GST reference and fragment maps for minibwa BAM analysis.")
    parser.add_argument("--reference", required=True)
    parser.add_argument("--metadata", required=True)
    parser.add_argument("--out-reference", required=True)
    parser.add_argument("--results-dir", default="analysis_minibwa_bam/results")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.out_reference), exist_ok=True)
    os.makedirs(args.results_dir, exist_ok=True)

    record = next(SeqIO.parse(args.reference, "fasta"))
    ref_id = record.id
    ref_seq = normalize(str(record.seq))
    with open(args.out_reference, "w") as handle:
        handle.write(f">{ref_id}\n")
        for i in range(0, len(ref_seq), 80):
            handle.write(ref_seq[i:i + 80] + "\n")

    codon_rows = []
    for codon_index, start0 in enumerate(range(0, len(ref_seq), 3), start=1):
        codon = ref_seq[start0:start0 + 3]
        if len(codon) < 3:
            aa = "partial"
        else:
            aa = CODON_TABLE.get(codon, "X")
        is_n_terminal_nls = 2 <= codon_index <= 8
        is_c_terminal_nls = 426 <= codon_index <= 432
        is_stop = aa == "*"
        if codon_index == 1:
            plot_pos = 1
        elif is_n_terminal_nls or is_c_terminal_nls or is_stop:
            plot_pos = ""
        else:
            plot_pos = codon_index - 7
        codon_rows.append({
            "reference_aa_position": codon_index,
            "plot_aa_position": plot_pos,
            "nt_start_1based": start0 + 1,
            "nt_end_1based": start0 + len(codon),
            "wt_codon": codon,
            "wt_aa": aa,
            "is_n_terminal_nls": is_n_terminal_nls,
            "is_c_terminal_nls": is_c_terminal_nls,
            "is_stop": is_stop,
        })

    fragment_rows = []
    primer_rows = []
    for row in csv.DictReader(open(args.metadata, newline="")):
        sample = row["SampleNumber"]
        pcr = normalize(row["PCR_product"])
        pcr_start = ref_seq.find(pcr)
        if pcr_start < 0:
            raise ValueError(f"{sample} PCR_product does not map exactly to reference")
        pcr_end = pcr_start + len(pcr)
        f_primer = normalize(row["F_primer"])
        r_primer = normalize(row["R_primer"])
        r_primer_rc = str(Seq(r_primer).reverse_complement()).upper()
        f_start = ref_seq.find(f_primer)
        r_start = ref_seq.find(r_primer_rc)
        fragment_rows.append({
            "SampleNumber": sample,
            "library_round": row["library_round"],
            "fragment_id": sample.split("-")[-1],
            "fragment_start0": pcr_start,
            "fragment_end0": pcr_end,
            "fragment_start_1based": pcr_start + 1,
            "fragment_end_1based": pcr_end,
            "fragment_length": len(pcr),
            "first_codon": pcr_start // 3 + 1,
            "last_codon": (pcr_end - 1) // 3 + 1,
        })
        primer_rows.append({
            "SampleNumber": sample,
            "F_primer_maps": f_start >= 0,
            "F_primer_start_1based": f_start + 1 if f_start >= 0 else "",
            "F_primer_matches_fragment_5prime": pcr.startswith(f_primer),
            "R_primer_rc_maps": r_start >= 0,
            "R_primer_rc_start_1based": r_start + 1 if r_start >= 0 else "",
            "R_primer_rc_matches_fragment_3prime": pcr.endswith(r_primer_rc),
        })

    write_csv(
        os.path.join(args.results_dir, "reference_summary.csv"),
        ["reference_id", "reference_length_bp", "codon_count", "length_mod_3", "terminal_codon"],
        [{
            "reference_id": ref_id,
            "reference_length_bp": len(ref_seq),
            "codon_count": len(codon_rows),
            "length_mod_3": len(ref_seq) % 3,
            "terminal_codon": ref_seq[-3:],
        }],
    )
    write_csv(os.path.join(args.results_dir, "reference_codon_table.csv"), list(codon_rows[0]), codon_rows)
    write_csv(os.path.join(args.results_dir, "fragment_map.csv"), list(fragment_rows[0]), fragment_rows)
    write_csv(os.path.join(args.results_dir, "primer_map_qc.csv"), list(primer_rows[0]), primer_rows)


if __name__ == "__main__":
    main()
