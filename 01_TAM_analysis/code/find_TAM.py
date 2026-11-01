import argparse
import csv
import logging
import gzip
import itertools
from typing import Dict

def count_subsequences(input_file: str, output_file: str, target_sequence: str, length: int) -> None:
    target_sequence = target_sequence.upper()  # Convert to uppercase
    if not target_sequence or len(target_sequence) < 6:
        raise ValueError("Target sequence must be at least 6 bases long.")
    if length < 1:
        raise ValueError("Length must be at least 1.")

    counts: Dict[str, int] = {}
    
    logging.basicConfig(level=logging.INFO)
    logging.info(f"Processing file: {input_file}")
    
    try:
        # Choose the appropriate opener for the file type
        open_func = gzip.open if input_file.endswith('.gz') else open
        
        with open_func(input_file, 'rt', encoding='utf-8') as f:
            for record in itertools.zip_longest(*[f] * 4):  # Each record contains four lines
                sequence = record[1].strip()  # The second line contains the sequence

                # Check whether the sequence contains valid characters
                if not all(c in 'ATCGN' for c in sequence):  # Allow only valid bases
                    logging.warning(f"Invalid characters found in sequence: {sequence}")
                    continue

                target_index = sequence.find(target_sequence)

                if target_index != -1 and target_index + len(target_sequence) + length <= len(sequence):
                    subsequence = sequence[target_index + len(target_sequence):target_index + len(target_sequence) + length]
                    
                    # Check whether the subsequence contains 'N'
                    if 'N' in subsequence:
                        logging.warning(f"Subsequence contains 'N', skipping: {subsequence}")
                        continue
                    
                    counts[subsequence] = counts.get(subsequence, 0) + 1

    except FileNotFoundError:
        raise FileNotFoundError(f"Input file {input_file} not found.")
    except Exception as e:
        logging.error(f"Error processing file: {e}")
        raise

    total_count = sum(counts.values())  # Count all retained subsequences
    
    # Write the results to a CSV file
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Subsequence', 'Count', 'Proportion'])  # Column names
        
        if counts:
            for seq, count in counts.items():
                proportion = count / total_count  # Calculate the proportion
                writer.writerow([seq, count, f"{proportion:.6f}"])  # Format to six decimal places
        else:
            logging.info("No valid subsequences found.")

    logging.info(f"Results written to: {output_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Count subsequences following a target sequence in a FASTA or FASTQ file.")
    parser.add_argument('--input', required=True, help='The input FASTA or FASTQ file.')
    parser.add_argument('--output', required=True, help='The output count file.')
    parser.add_argument('--target_seq', required=True, help='The target sequence to search for.')
    parser.add_argument('--length', type=int, default=6, help='Length of bases to find after the target sequence (default: 6).')  # Subsequence length

    args = parser.parse_args()

    count_subsequences(args.input, args.output, args.target_seq, args.length)
