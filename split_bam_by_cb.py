#!/usr/bin/env python3
import os
import argparse
import pysam

def split_bam_by_cb(input_bam, output_dir):
    """
    Splits a BAM file based on the 'CB' tag and writes separate BAM files for each barcode.

    Parameters:
      input_bam: Path to the input BAM file.
      output_dir: Directory where the split BAM files will be stored.
    """
    # Open the input BAM file
    bam_in = pysam.AlignmentFile(input_bam, "rb")
    
    # Dictionary to store output BAM file handles for each unique CB tag
    output_handles = {}
    
    # Iterate through each alignment in the BAM file
    for read in bam_in.fetch(until_eof=True):
        try:
            cb = read.get_tag("CB")
        except KeyError:
            # Skip records that do not have a 'CB' tag
            continue

        # Rename barcode '-' to 'noCB' to avoid ambiguity in bash
        if cb == "-":
            cb = "noCB"
        
        # If this barcode is not yet encountered, create a new output BAM file
        if cb not in output_handles:
            output_path = os.path.join(output_dir, f"{cb}.bam")
            output_handles[cb] = pysam.AlignmentFile(output_path, "wb", template=bam_in)
        output_handles[cb].write(read)
    
    # Close the input BAM file
    bam_in.close()
    
    # Close all output BAM file handles
    for handle in output_handles.values():
        handle.close()

def main():
    parser = argparse.ArgumentParser(
        description="Split a BAM file based on the 'CB' tag using pysam. Supports specifying input file and output directory."
    )
    parser.add_argument("-i", "--input", required=True, help="Path to the input BAM file")
    parser.add_argument("-o", "--output_dir", required=True, help="Path to the output directory")
    
    args = parser.parse_args()
    
    # Check if the input file exists
    if not os.path.exists(args.input):
        parser.error(f"Input file {args.input} does not exist.")
    
    # Create the output directory if it does not exist
    if not os.path.isdir(args.output_dir):
        os.makedirs(args.output_dir, exist_ok=True)
    
    split_bam_by_cb(args.input, args.output_dir)
    
if __name__ == '__main__':
    main()
