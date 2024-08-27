#!/usr/bin/env python2
# -*- coding: utf-8 -*-

"""
Author: Hui Zhou
Version: 1.1
Date: 2024-08-21
Description: This script calculates and plots base frequencies in a FASTQ file or from an existing CSV table.
             Outputs include CSV table, static image (PNG by default), and optional dynamic interactive HTML chart.
Usage:
    module load python (python/3.12.1)
    - To analyze a FASTQ file and generate output files:
        python script.py -i example.fastq -o output_prefix
    - Specify bases at specific positions, generate dynamic interactive HTML, and add base percentage to PNG:
        python script.py -i example.fastq -o output_prefix -s 8 16 --interactive --add_percentage -sp special_output_prefix
    - Generate plots directly from an existing CSV table:
        python script.py --table_input existing_frequencies.csv -o output_prefix --interactive
"""

import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use the Agg backend for non-GUI environments
import matplotlib.pyplot as plt
import plotly.graph_objs as go
from Bio import SeqIO
import argparse
import gzip
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
    Handles both individual positions and ranges of positions.
    """
    if len(positions) == 2:
        start_pos, end_pos = positions
        positions_data = df.iloc[start_pos:end_pos+1]
    else:
        positions_data = df.iloc[positions]
    return positions_data

def add_percentage_labels(ax, df, colors):
    """
    Adds percentage labels below each base position and in the same color as the line.
    """
    for i in range(len(df)):
        if (i + 1) % 5 == 0:  # Optional: Only add labels for positions that are multiples of 5
            labels = []
            for base, color in zip(df.columns, colors):
                percentage = df.iloc[i][base] * 100
                if i == 0:
                    label = '{}{:0.0f}%'.format(base, percentage)
                else:
                    label = '{:0.0f}%'.format(percentage)
                ax.text(i + 1, -0.1 - 0.05 * df.columns.get_loc(base), label, ha='center', va='top', 
                        fontsize=8, color=color, transform=ax.get_xaxis_transform())

def plot_base_frequencies_static(df, output_plot, file_name, n_bases, add_percentage):
    """
    Plots the base frequencies as a static image.
    """
    plt.figure(figsize=(12, 8))
    ax = plt.gca()
    colors = ['blue', 'orange', 'green', 'red', 'purple']  # Specify the colors for A, T, G, C, N
    for base, color in zip(df.columns, colors):
        ax.plot(df.index + 1, df[base], label=base, color=color)
    
    if add_percentage:
        add_percentage_labels(ax, df, colors)
    
    plt.title('Base Frequencies in {} (First {} Bases)'.format(file_name, n_bases))
    plt.xlabel('Base Position')
    plt.ylabel('Frequency')
    plt.legend(loc='upper right')
    plt.grid(True)

    # Set x-axis to show integer values only
    plt.xticks(range(1, len(df) + 1, 5))

    # Adjust the bottom margin to fit the labels if percentage labels are added
    if add_percentage:
        plt.subplots_adjust(bottom=0.3)

    plt.savefig(output_plot)
    plt.close()

def plot_specific_positions_frequencies(df, output_plot, file_name, start_pos, end_pos):
    """
    Plots the base frequencies for specific positions as a static image.
    """
    plt.figure(figsize=(12, 8))
    ax = plt.gca()
    colors = ['blue', 'orange', 'green', 'red', 'purple']  # Specify the colors for A, T, G, C, N
    for base, color in zip(df.columns, colors):
        ax.plot(df.index + 1, df[base], label=base, color=color)
    
    plt.title('Base Frequencies at Specific Positions from {} to {} in {}'.format(start_pos, end_pos, file_name))
    plt.xlabel('Base Position')
    plt.ylabel('Frequency')
    plt.legend(loc='upper right')
    plt.grid(True)

    plt.xticks(range(start_pos, end_pos + 1, 1))

    plt.savefig(output_plot)
    plt.close()

def plot_base_frequencies_interactive(df, output_plot, file_name, n_bases):
    """
    Plots the base frequencies as an interactive Plotly graph.
    """
    fig = go.Figure()
    colors = ['blue', 'orange', 'green', 'red', 'purple']  # Specify the colors for A, T, G, C, N
    
    for base, color in zip(df.columns, colors):
        fig.add_trace(go.Scatter(
            x=df.index + 1,
            y=df[base],
            mode='lines+markers',
            name=base,
            line=dict(color=color),
            text=['Position {}: {}: {:.2%}'.format(i + 1, base, v) for i, v in enumerate(df[base])],  # Add position and percentage as hover text
            hoverinfo='text'
        ))

    fig.update_layout(
        title='Base Frequencies in {} (First {} Bases)'.format(file_name, n_bases),
        xaxis_title='Base Position',
        yaxis_title='Frequency',
        legend_title='Base',
        template='plotly_white'
    )
    
    fig.write_html(output_plot)

def main(input_file, table_input, output_prefix, n_bases, specific_positions, interactive, add_percentage, specific_output_prefix):
    # Determine file name and output paths
    file_name = os.path.basename(input_file.split(".")[0]) if input_file else os.path.basename(table_input.split(".")[0])
    output_table = "{}.csv".format(output_prefix)
    output_plot = "{}.png".format(output_prefix)
    output_html = "{}.html".format(output_prefix)

    # Parse input file or load existing table
    if table_input:
        base_frequencies = pd.read_csv(table_input, index_col=0)
        if n_bases < len(base_frequencies):
            base_frequencies = base_frequencies.iloc[:n_bases]
    else:
        base_frequencies = parse_fastq(input_file, n_bases)
        base_frequencies.to_csv(output_table, index_label='Base Position')

    # Plot specific positions if requested
    if specific_positions:
        specific_positions = [pos - 1 for pos in specific_positions]  # Convert to zero-based indexing
        positions_data = extract_specific_positions(base_frequencies, specific_positions)
        start_pos, end_pos = specific_positions[0] + 1, specific_positions[-1] + 1  # Convert to 1-based indexing for display
        specific_output_plot = "{}.png".format(specific_output_prefix)
        plot_specific_positions_frequencies(positions_data, specific_output_plot, file_name, start_pos, end_pos)

    # Generate static plot
    plot_base_frequencies_static(base_frequencies, output_plot, file_name, n_bases, add_percentage)

    # Generate interactive plot if requested
    if interactive:
        plot_base_frequencies_interactive(base_frequencies, output_html, file_name, n_bases)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate and plot base frequencies in a FASTQ file or from an existing CSV table.')
    parser.add_argument('-i', '--input_file', type=str, help='Path to the FASTQ file, support .gz')
    parser.add_argument('-a', '--table_input', type=str, help='Path to an existing base frequencies table to use as input instead of a FASTQ file.')
    parser.add_argument('-o', '--output_prefix', type=str, required=True, help='Output file name prefix.')
    parser.add_argument('-n', '--num_bases', type=int, default=85, help='Number of bases to analyze (default: 85).')
    parser.add_argument('-s', '--specific_positions', type=int, nargs='+', metavar='POS', help='Specific base positions or range to analyze (e.g., -s 10 15)')
    parser.add_argument('-sp', '--specific_output_prefix', type=str, help='Output file name prefix for specific positions plot.')
    parser.add_argument('--interactive', action='store_true', help='Generate an interactive HTML plot.')
    parser.add_argument('--add_percentage', action='store_true', help='Add percentage labels to the static image.')


    args = parser.parse_args()
    
    main(args.input_file, args.table_input, args.output_prefix, args.num_bases, args.specific_positions, args.interactive, args.add_percentage, args.specific_output_prefix)
