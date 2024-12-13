#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#module load python/3.12.3
"""
Author: Hui Zhou
Version: 0.1.0
Date: 2024-12-13
"""

import argparse
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
from itertools import combinations
import logging
import numpy as np
import matplotlib.pyplot as plt

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

def hamming_distance(seq1, seq2):
    """Calculate the Hamming distance between two sequences."""
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of equal length to calculate Hamming distance.")
    return sum(el1 != el2 for el1, el2 in zip(seq1, seq2))

def sliding_window(seq, window_size):
    """Generate sliding windows of a given size from a sequence."""
    return [seq[i:i+window_size] for i in range(len(seq) - window_size + 1)]

def calculate_hamming_matrix(barcodes):
    """Calculate the Hamming distance matrix for a list of barcodes."""
    num_barcodes = len(barcodes)
    matrix = np.zeros((num_barcodes, num_barcodes), dtype=int)
    for i, record1 in enumerate(barcodes):
        for j, record2 in enumerate(barcodes):
            if i < j:
                matrix[i, j] = hamming_distance(str(record1.seq), str(record2.seq))
                matrix[j, i] = matrix[i, j]  # Symmetry
    return matrix

def save_hamming_matrix(matrix, barcodes, output_file):
    """Save the Hamming distance matrix to a text file."""
    with open(output_file, "w") as f:
        f.write("\t" + "\t".join(record.id for record in barcodes) + "\n")
        for i, record in enumerate(barcodes):
            f.write(record.id + "\t" + "\t".join(map(str, matrix[i])) + "\n")

def plot_hamming_matrix(matrix, barcodes, output_image):
    """Plot the Hamming distance matrix as a heatmap."""
    plt.figure(figsize=(10, 8))
    plt.imshow(matrix, cmap="viridis", interpolation="nearest")
    plt.colorbar(label="Hamming Distance")
    plt.xticks(ticks=range(len(barcodes)), labels=[record.id for record in barcodes], rotation=90)
    plt.yticks(ticks=range(len(barcodes)), labels=[record.id for record in barcodes])
    plt.title("Hamming Distance Matrix")
    plt.tight_layout()
    plt.savefig(output_image)
    plt.close()

def filter_barcodes(input_fasta, exclude_fasta, hamming_threshold, gc_min, gc_max, output_filtered, output_info, hamming_matrix_file, hamming_plot_file, debug):
    """Filter barcodes based on Hamming distance, GC content, and exclusions."""
    if debug:
        logging.getLogger().setLevel(logging.DEBUG)

    # Read barcodes from input fasta
    barcodes = list(SeqIO.parse(input_fasta, "fasta"))
    logging.debug(f"Loaded {len(barcodes)} barcodes from input file.")

    # Read exclusion sequences from exclude fasta
    exclusions = list(SeqIO.parse(exclude_fasta, "fasta")) if exclude_fasta else []
    logging.debug(f"Loaded {len(exclusions)} exclusion sequences from exclude file.")

    # Prepare results
    filtered_barcodes = []
    info_lines = ["ID\tSequence\tGC_Content\tExcluded_Match\tMatching_Exclusion_Sequence"]

    for record in barcodes:
        seq = str(record.seq)
        gc_content = gc_fraction(seq) * 100

        # Check GC content
        if gc_content < gc_min or gc_content > gc_max:
            logging.debug(f"Barcode {record.id} excluded due to GC content: {gc_content:.2f}%")
            info_lines.append(f"{record.id}\t{seq}\t{gc_content:.2f}\tGC_Content_Out_Range\tNA")
            continue

        # Check Hamming distance with other barcodes
        hamming_pass = True
        for other_record in filtered_barcodes:
            if hamming_distance(seq, str(other_record.seq)) < hamming_threshold:
                logging.debug(f"Barcode {record.id} excluded due to low Hamming distance with {other_record.id}.")
                info_lines.append(f"{record.id}\t{seq}\t{gc_content:.2f}\tLow_Hamming_Distance\tNA")
                hamming_pass = False
                break

        if not hamming_pass:
            continue

        # Check against exclusion sequences using sliding window
        exclude_match = False
        matching_exclusion = ""
        for exclude_record in exclusions:
            exclude_seq = str(exclude_record.seq)
            for window in sliding_window(exclude_seq, len(seq)):
                if seq == window or seq == str(window[::-1].translate(str.maketrans("ATCG", "TAGC"))):
                    logging.debug(f"Barcode {record.id} excluded due to match with exclusion sequence {exclude_record.id}.")
                    matching_exclusion = exclude_seq
                    info_lines.append(f"{record.id}\t{seq}\t{gc_content:.2f}\tExcluded_Match\t{matching_exclusion}")
                    exclude_match = True
                    break
            if exclude_match:
                break

        if exclude_match:
            continue

        # If passed all filters, add to filtered barcodes
        filtered_barcodes.append(record)
        logging.debug(f"Barcode {record.id} passed all filters.")
        info_lines.append(f"{record.id}\t{seq}\t{gc_content:.2f}\tPassed\tNA")

    # Write filtered barcodes to output
    with open(output_filtered, "w") as filtered_handle:
        SeqIO.write(filtered_barcodes, filtered_handle, "fasta")
    logging.debug(f"Filtered barcodes written to {output_filtered}.")

    # Write detailed info to output
    with open(output_info, "w") as info_handle:
        info_handle.write("\n".join(info_lines))
    logging.debug(f"Filtering details written to {output_info}.")

    # Calculate Hamming distance matrix
    hamming_matrix = calculate_hamming_matrix(barcodes)
    save_hamming_matrix(hamming_matrix, barcodes, hamming_matrix_file)
    logging.debug(f"Hamming matrix written to {hamming_matrix_file}.")

    # Plot Hamming distance matrix
    plot_hamming_matrix(hamming_matrix, barcodes, hamming_plot_file)
    logging.debug(f"Hamming matrix plot saved to {hamming_plot_file}.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter barcode sequences based on Hamming distance, GC content, and exclusions.")
    parser.add_argument("-i", "--input", required=True, help="Input fasta file containing barcodes.")
    parser.add_argument("-e", "--exclude", required=False, help="Fasta file containing sequences to exclude.")
    parser.add_argument("-ht", "--hamming", type=int, default=3, help="Minimum Hamming distance threshold (default: 3).")
    parser.add_argument("-gmin", "--gc_min", type=float, default=40.0, help="Minimum GC content percentage (default: 40.0).")
    parser.add_argument("-gmax", "--gc_max", type=float, default=60.0, help="Maximum GC content percentage (default: 60.0).")
    parser.add_argument("-of", "--output_filtered", required=True, help="Output fasta file for filtered barcodes.")
    parser.add_argument("-oi", "--output_info", required=True, help="Output text file for detailed filtering information.")
    parser.add_argument("-hmf", "--hamming_matrix_file", required=True, help="Output text file for Hamming distance matrix.")
    parser.add_argument("-hpf", "--hamming_plot_file", required=True, help="Output image file for Hamming distance heatmap.")
    parser.add_argument("-d", "--debug", action="store_true", help="Enable debug logging.")

    args = parser.parse_args()

    filter_barcodes(
        input_fasta=args.input,
        exclude_fasta=args.exclude,
        hamming_threshold=args.hamming,
        gc_min=args.gc_min,
        gc_max=args.gc_max,
        output_filtered=args.output_filtered,
        output_info=args.output_info,
        hamming_matrix_file=args.hamming_matrix_file,
        hamming_plot_file=args.hamming_plot_file,
        debug=args.debug
    )




