#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#module load python/3.12.3

"""
Author: Hui Zhou
Version: 0.2.0
Date: 2025-Jul-7

Description:
    Remove soft-clipped bases from both ends of reads in a BAM file.
    If the 5' soft clip exactly matches a user-specified pattern (supports IUPAC codes),
    also output the trimmed read to a separate BAM file.
    Additionally, output statistics for all observed 5' and 3' softclip contents.
"""

import pysam
import argparse
from collections import Counter

# IUPAC codes for DNA bases
IUPAC = {
    'A': 'A',
    'C': 'C',
    'G': 'G',
    'T': 'T',
    'R': 'AG',
    'Y': 'CT',
    'S': 'GC',
    'W': 'AT',
    'K': 'GT',
    'M': 'AC',
    'B': 'CGT',
    'D': 'AGT',
    'H': 'ACT',
    'V': 'ACG',
    'N': 'ACGT',
}

def seq_match_pattern(seq, pattern):
    """Compare sequence to pattern with IUPAC code support."""
    if len(seq) != len(pattern):
        return False
    for s, p in zip(seq, pattern):
        if s.upper() not in IUPAC.get(p.upper(), p.upper()):
            return False
    return True

def remove_softclip_and_tag(read, pattern):
    """
    Remove soft clips at both ends of the read.
    If the original 5' softclip region matches the pattern, return flag True.
    Return:
        trimmed_read: pysam.AlignedSegment (modified in-place)
        matched_flag: bool (True if matched)
        left_seq: str (original 5' softclip bases, '' if none)
        right_seq: str (original 3' softclip bases, '' if none)
    """
    cigar = read.cigartuples
    seq = read.query_sequence
    qual = read.query_qualities

    left_clip = 0
    right_clip = 0
    left_seq = ""
    right_seq = ""

    # Check 5' soft clip
    if cigar and cigar[0][0] == 4:  # 4 = S
        left_clip = cigar[0][1]
        left_seq = seq[:left_clip]
    # Check 3' soft clip
    if cigar and cigar[-1][0] == 4:
        right_clip = cigar[-1][1]
        right_seq = seq[-right_clip:]
    # Remove soft clips
    new_seq = seq[left_clip:len(seq)-right_clip if right_clip != 0 else None]
    new_qual = qual[left_clip:len(qual)-right_clip if right_clip != 0 else None]
    # Build new CIGAR without S at ends
    new_cigar = []
    for idx, (op, length) in enumerate(cigar):
        if (idx == 0 or idx == len(cigar)-1) and op == 4:
            continue
        new_cigar.append((op, length))
    # Modify read in-place
    read.query_sequence = new_seq
    read.query_qualities = new_qual
    read.cigartuples = new_cigar
    # Only match if left_clip == len(pattern) and matches pattern
    matched_flag = (left_clip == len(pattern)) and seq_match_pattern(left_seq, pattern)
    return read, matched_flag, left_seq, right_seq

def write_softclip_stats(left_counter, right_counter, total_reads, out_path):
    """Write softclip stats (counts and frequencies) to a text file."""
    with open(out_path, 'w') as out:
        out.write("5' softclip statistics:\n")
        out.write(f"{'Softclip':<20}\t{'Count':<10}\t{'Frequency'}\n")
        for seq, count in left_counter.most_common():
            freq = count / total_reads
            out.write(f"{seq:<20}\t{count:<10}\t{freq:.6f}\n")
        out.write("\n3' softclip statistics:\n")
        out.write(f"{'Softclip':<20}\t{'Count':<10}\t{'Frequency'}\n")
        for seq, count in right_counter.most_common():
            freq = count / total_reads
            out.write(f"{seq:<20}\t{count:<10}\t{freq:.6f}\n")

def main():
    parser = argparse.ArgumentParser(
        description="Remove softclipped bases from BAM. "
                    "If 5' softclip exactly matches pattern (IUPAC allowed), output to separate BAM. "
                    "Additionally, collect statistics for all 5' and 3' softclip sequences."
    )
    parser.add_argument('-i', '--input', required=True, help='Input BAM file')
    parser.add_argument('-o', '--output', required=True, help='Output BAM file (softclip removed)')
    parser.add_argument('--pattern_bam', required=True, help='Output BAM for reads with 5\' softclip matching pattern')
    parser.add_argument('--trim_pattern', required=True, help='5\' pattern to match (e.g. GGGAC, WWGGG, etc.)')
    parser.add_argument('--softclip_stats', required=True, help='Output txt file for softclip statistics')
    args = parser.parse_args()

    pattern = args.trim_pattern.upper()

    left_counter = Counter()
    right_counter = Counter()
    total_reads = 0

    with pysam.AlignmentFile(args.input, "rb") as infile, \
         pysam.AlignmentFile(args.output, "wb", header=infile.header) as out_bam, \
         pysam.AlignmentFile(args.pattern_bam, "wb", header=infile.header) as pattern_bam:
        for read in infile:
            total_reads += 1
            if read.is_unmapped or read.cigartuples is None:
                out_bam.write(read)
                continue
            trimmed_read, matched, left_seq, right_seq = remove_softclip_and_tag(read, pattern)
            out_bam.write(trimmed_read)
            if matched:
                pattern_bam.write(trimmed_read)
            # Count softclips (use '' if no softclip, for completeness)
            left_counter[left_seq] += 1
            right_counter[right_seq] += 1

    write_softclip_stats(left_counter, right_counter, total_reads, args.softclip_stats)

if __name__ == "__main__":
    main()
