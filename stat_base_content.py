#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Author: Hui Zhou
Version: 1.0
Date: 2024-08-21
Description: This script calculates and plots base frequencies in a FASTQ file. It can also specify base positions of interest. It also supports plotting from an existing frequency table.
Usage:
    - To analyze a FASTQ file and generate frequency plots:
        python script.py -i example.fastq -n 85 -t output_table.csv -p plot.png
    - Specify bases at specific positions simultaneously, for example, from 8 to 16:
        python script.py -i example.fastq -n 85 -t output_table.csv -p plot.png -s 8 16 -sp 8_16_base.png
    - To generate plots directly from an existing frequency table:
        python script.py --table_input existing_frequencies.csv -p plot.png -s 10 15 -sp custom_name.png
"""

import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use the Agg backend for non-GUI environments
import matplotlib.pyplot as plt
from Bio import SeqIO
import numpy as np
import argparse
import gzip
import sys
import os

def parse_fastq(file_path, n_bases):
    base_counts = {base: [0] * n_bases for base in 'ATCGN'}
    total_counts = [0] * n_bases

    # Detect if the file is gzipped
    if file_path.endswith('.gz'):
        handle = gzip.open(file_path, 'rt')  # Open as text
    else:
        handle = open(file_path, 'r')

    with handle:
        for record in SeqIO.parse(handle, 'fastq'):
            sequence = str(record.seq)
            for i in range(min(n_bases, len(sequence))):
                base = sequence[i]
                if base in base_counts:
                    base_counts[base][i] += 1
                    total_counts[i] += 1

    # Calculate frequencies with explicit float conversion
    base_freqs = {
        base: [count / float(total) if total > 0 else 0 for count, total in zip(counts, total_counts)]
        for base, counts in base_counts.items()
    }

    return pd.DataFrame(base_freqs)

def extract_specific_positions(df, positions):
    """
    Extracts the frequencies of specific base positions from the DataFrame.
    """
    positions_data = df.iloc[positions]
    return positions_data

def plot_base_frequencies(df, output_plot, file_name, n_bases):
    """
    Plots the base frequencies and ensures the x-axis uses integer values.
    The plot title includes the input file name and the number of bases analyzed.
    """
    plt.figure(figsize=(10, 6))
    for base in df.columns:
        plt.plot(df.index + 1, df[base], label=base)
    
    plt.title('Base Frequencies in {} (First {} Bases)'.format(file_name, n_bases))
    plt.xlabel('Base Position')
    plt.ylabel('Frequency')
    plt.legend()
    plt.grid(True)

    # Set x-axis to show integer values only
    plt.xticks(range(0, len(df) + 1, 5))

    plt.savefig(output_plot)
    plt.close()

def plot_specific_positions_frequencies(df, positions, output_plot, file_name):
    """
    Plots the base frequencies at specific positions.
    """
    plt.figure(figsize=(10, 6))
    for base in df.columns:
        plt.plot([pos + 1 for pos in positions], df.loc[positions, base], marker='o', label=base)
    
    plt.title('Base Frequencies at Specific Positions in {}'.format(file_name))
    plt.xlabel('Base Position')
    plt.ylabel('Frequency')
    plt.legend()
    plt.grid(True)

    plt.savefig(output_plot)
    plt.close()

def main(input_file, n_bases, output_table, output_plot, specific_positions, specific_plot_name, table_input):
    
    if table_input:
        base_frequencies = pd.read_csv(table_input, index_col=0)
        file_name = os.path.basename(table_input.split(".")[0])
        
        if n_bases < len(base_frequencies):
            base_frequencies = base_frequencies.iloc[:n_bases]
    else:
        # Extract file name without path for the plot title
        file_name = os.path.basename(input_file.split(".")[0])

        # Parse the FASTQ file and calculate base frequencies
        base_frequencies = parse_fastq(input_file, n_bases)
        
        # Save the frequencies table to a file
        base_frequencies.to_csv(output_table, index_label='Base Position')

    # Extract and plot specific positions if provided
    if specific_positions:
        specific_positions = [pos - 1 for pos in specific_positions]  # Convert to zero-based indexing
        positions_data = extract_specific_positions(base_frequencies, specific_positions)
        print("Frequencies at specified positions:")
        print(positions_data)

        # Use the custom name if provided, otherwise generate a default name
        if specific_plot_name:
            specific_plot_file = specific_plot_name
               
        plot_specific_positions_frequencies(base_frequencies, specific_positions, specific_plot_file, file_name)

    # Plot the base frequencies and save the plot
    plot_base_frequencies(base_frequencies, output_plot, file_name, n_bases)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate and plot base frequencies in a FASTQ file or from an existing frequency table.')
    parser.add_argument('-i', '--input_file', type=str, help='Path to the FASTQ file, support .gz')
    parser.add_argument('-n', '--num_bases', type=int, default=85, help='Number of bases to analyze (default: 85).')
    parser.add_argument('-t', '--output_table', type=str, required=False, help='Output table file path.')
    parser.add_argument('-p', '--output_plot', type=str, required=True, help='Output plot file path.')
    parser.add_argument('-s', '--specific_positions', type=int, nargs='+', metavar='POS', help='Specific base positions to analyze (e.g., -s 10 15)')
    parser.add_argument('-sp','--specific_plot_name', type=str, help='Custom name for the specific positions plot output.')
    parser.add_argument('-a','--table_input', type=str, help='Path to an existing base frequencies table to use as input instead of a FASTQ file.')

    args = parser.parse_args()
    
    main(args.input_file, args.num_bases, args.output_table, args.output_plot, args.specific_positions, args.specific_plot_name, args.table_input)