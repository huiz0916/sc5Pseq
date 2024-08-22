#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Author: Hui Zhou
Version: 1.0
Date: 2024-08-19
"""

import os
import sys
import pandas as pd

def parse_final_out(file_path):
    """
    Parses the .final.out file and extracts the relevant information.
    """
    with open(file_path, 'r') as f:
        lines = f.readlines()

    # Initialize a list to store the relevant data
    parsed_data = []
    keys_order = []

    # Filter out lines that are not needed and clean the remaining lines
    for line in lines:
        line = line.strip()  # Remove leading/trailing whitespaces
        if line.startswith("UNIQUE READS:") or line.startswith("MULTI-MAPPING READS:") or \
           line.startswith("UNMAPPED READS:") or line.startswith("CHIMERIC READS:"):
            continue  # Skip the unwanted lines

        if "|" in line:
            key, value = line.split("|")
            key = key.strip()  # Remove any leading/trailing whitespaces
            value = value.strip()  # Remove any leading/trailing whitespaces
            parsed_data.append((key, value))
            keys_order.append(key)

    return dict(parsed_data), keys_order

def process_directory(directory_path):
    """
    Processes all .final.out files in the given directory and returns a combined DataFrame.
    """
    combined_data = {}
    columns_order = None

    for filename in os.listdir(directory_path):
        if filename.endswith(".final.out"):
            file_path = os.path.join(directory_path, filename)
            file_data, keys_order = parse_final_out(file_path)
            # Use the filename (without extension) as the column name
            column_name = os.path.splitext(filename)[0]
            combined_data[column_name] = file_data
            # Store the keys order based on the first file processed
            if columns_order is None:
                columns_order = keys_order

    # Convert the dictionary of data into a DataFrame
    df = pd.DataFrame(combined_data)
    df.index.name = 'Metrics'

    # Reorder the DataFrame's index to match the original file's order
    df = df.reindex(columns_order)

    return df

def main(directory_path, output_name):
    df = process_directory(directory_path)
    # Save the combined DataFrame to a CSV file
    df.to_csv(output_name, index=True)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py DirectoryPath OutputName")
    else:
        directory_path = sys.argv[1]
        output_name = sys.argv[2]
        main(directory_path, output_name)
