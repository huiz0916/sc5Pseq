#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Author: Hui Zhou
Version: 1.0
Date: 2024-12-06
"""

"""
Easily plot the format mapping stat results from STAR.

Why? Sometimes I have many samples to check the mapping rate, and it's not easy to visualize the sample mapping variation.

Usage:
python star_plot.py -i [file] -f [specified features, e.g., umr:Uniquely mapped reads] -o [Output file name prefix] --interactive -c [color]
"""
import os
import sys
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import plotly.express as px


FEATURE_ABBREVIATIONS = {
    "Number of input reads": "input",
    "Uniquely mapped reads number": "umr_num",
    "Uniquely mapped reads %": "umr_pct",
    "Mismatch rate per base, %": "mismatch_rate",
    "Average mapped length": "avg_len",
    "Number of reads mapped to multiple loci": "multi_num",
    "% of reads mapped to multiple loci": "multi_pct",
    "Number of reads unmapped: too short": "unmapped_short_num",
    "% of reads unmapped: too short": "unmapped_short_pct"
}

FEATURES_REVERSE_MAP = {v: k for k, v in FEATURE_ABBREVIATIONS.items()}

# Science classic colors
SCIENCE_COLORS = {
    'blue': '#377eb8',
    'orange': '#ff7f00',
    'green': '#4daf4a',
    'red': '#e41a1c',
    'purple': '#984ea3'
}
DEFAULT_COLOR = 'blue'


def plot_features(df, features, output_prefix, interactive, color):
    """
    Plot the specified features for the samples in the dataframe.
    
    Args:
        df (pd.DataFrame): The dataframe containing the data.
        features (list): List of features to plot.
        output_prefix (str): Prefix for the output file names.
        interactive (bool): Whether to generate interactive plots using Plotly.
        color (str): The color to use for the plots.
    """
    # Iterate over each feature to plot
    for feature in features:
        # Map abbreviated feature to full feature name if necessary
        full_feature = FEATURES_REVERSE_MAP.get(feature, feature)
        if full_feature not in df.index:
            print(f"Feature '{feature}' not found in the data. Skipping.")
            continue

        # Extract data for the current feature and handle percentage values
        data = df.loc[full_feature].str.replace('%', '').astype(float).dropna()
        if data.empty:
            print(f"Feature '{feature}' has no valid numeric data to plot. Skipping.")
            continue
        
        plt.figure(figsize=(10, 6))
        plt.title(f"Comparison of {full_feature}")
        plt.ylabel(full_feature)
        plt.xlabel("Sample Name")
        data.plot(kind='bar', color=color, label=full_feature)
        plt.xticks(rotation=45, ha='right')
        plt.legend()
        plt.tight_layout()
        # Save static plot
        plt.savefig(f"{output_prefix}_{feature}.png")
        plt.close()

        # Create interactive plot if requested
        if interactive:
            fig = px.bar(data.reset_index(), x='index', y=full_feature, title=f"Comparison of {full_feature}",
                         labels={'index': 'Sample Name', 'y': full_feature}, color_discrete_sequence=[color])
            fig.write_html(f"{output_prefix}_{feature}_interactive.html")


def main():
    """
    Main function to parse arguments and call the plotting function.
    """
    # Setup argument parser
    parser = argparse.ArgumentParser(description="Plot features from formatted STAR log data.")
    parser.add_argument("-i", "--input", required=True, help="Input CSV file containing formatted STAR log data.")
    parser.add_argument("-f", "--features", nargs='+', default=['umr_pct'],
                        help="List of features to compare (default is 'umr_pct'). Available features: input, umr_num, umr_pct, mismatch_rate, avg_len, multi_num, multi_pct, unmapped_short_num, unmapped_short_pct.")
    parser.add_argument("-c", "--color", choices=list(SCIENCE_COLORS.keys()), default=DEFAULT_COLOR,
                        help="Specify the color for the plots (default is 'blue'). Available colors: blue, orange, green, red, purple.")
    parser.add_argument("-o", "--output", required=True, help="Output file name prefix for the plots.")
    parser.add_argument("--interactive", action='store_true', help="Generate interactive plots using Plotly.")
    parser.add_argument("--sort-samples", action='store_true', help="Sort sample names in alphanumeric order.")


    args = parser.parse_args()

    # Load the data
    df = pd.read_csv(args.input, index_col='Metrics')

    # Remove 'Log.final' from sample names
    df.columns = df.columns.str.replace(r'Log\.final', '', regex=True)

    # Sort sample names if requested
    if args.sort_samples:
        df = df.reindex(sorted(df.columns), axis=1)

    # If no features specified, use the default feature
    if args.features is None:
        features = ['umr_pct']
    else:
        features = args.features

    # Plot the requested features
    plot_features(df, features, args.output, args.interactive, SCIENCE_COLORS[args.color])


if __name__ == "__main__":
    main()
