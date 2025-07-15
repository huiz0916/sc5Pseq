#!/usr/bin/env python3
# -*- coding: utf-8 -*-

#module load python/3.12.3
"""
Author: Hui Zhou
Version: 0.2.2
Date: 2025-Jul-11

Description:
    Pipeline to identify, visualize, and filter strand invasion/missing primer pairing artifacts in BAM files,
    based on hamming distance (excluding the 3' GGG/CCC of TSO) and upstream sequence G/C count.
    Includes:
        - stat:   Output TSVs for TSO-hamming+G, rcTSO-hamming+C, logo table (all reads), and grouped percentage curves
        - plot:   Visualization (heatmap, logo, or grouped curve) from stat outputs, supports flexible filtering
        - filtering: Filter BAM files by user-specified thresholds (strand_invasion/min_G, missing_pairing/min_C)
"""

import argparse
import pysam
from pyfaidx import Fasta
import pandas as pd
import matplotlib.pyplot as plt
import logomaker
from collections import Counter
import sys

IUPAC_CODES = {
    'A': {'A'}, 'C': {'C'}, 'G': {'G'}, 'T': {'T'},
    'R': {'A', 'G'}, 'Y': {'C', 'T'}, 'S': {'G', 'C'}, 'W': {'A', 'T'},
    'K': {'G', 'T'}, 'M': {'A', 'C'}, 'B': {'C', 'G', 'T'}, 'D': {'A', 'G', 'T'},
    'H': {'A', 'C', 'T'}, 'V': {'A', 'C', 'G'}, 'N': {'A', 'C', 'G', 'T'}
}

def iupac_match(base1, base2):
    """Return True if two bases match under IUPAC codes."""
    set1 = IUPAC_CODES.get(base1.upper(), set())
    set2 = IUPAC_CODES.get(base2.upper(), set())
    return bool(set1 & set2)

def hamming_distance_iupac(seq1, seq2):
    """Calculate Hamming distance with IUPAC codes."""
    assert len(seq1) == len(seq2)
    return sum(not iupac_match(b1, b2) for b1, b2 in zip(seq1, seq2))

def count_last3g(seq):
    """Count number of G in last 3 bases of seq."""
    return seq[-3:].upper().count('G')

def count_last3c(seq):
    """Count number of C in last 3 bases of seq."""
    return seq[-3:].upper().count('C')

def reverse_complement(seq):
    """Return reverse complement of a DNA sequence (IUPAC codes supported)."""
    complement = str.maketrans('ATCGNatcgnRYSWKMBDHV', 'TAGCNtagcnRYSWMKVHDB')
    return seq.translate(complement)[::-1]

def fetch_upstream_seq(fasta, chrom, pos, strand, length):
    """Fetch upstream genomic sequence of given length from read's 5' mapping site."""
    if strand == "+":
        start = max(pos - length, 0)
        end = pos
        seq = fasta[chrom][start:end].seq
    else:
        start = pos + 1
        end = pos + 1 + length
        seq = reverse_complement(fasta[chrom][start:end].seq)
    return seq

def get_core_tso(tsoseq, g3_len=3):
    """
    Extract core TSO (without the trailing GGG) and the G part length.
    Example: ATATGGG --> core=ATAT, GGG
    """
    core = tsoseq[:-g3_len]
    g3 = tsoseq[-g3_len:]
    return core, g3

def stat_mode(bam_path, fasta_path, tso_seq, g_out, c_out, logo_table, g_curve_out=None, c_curve_out=None):
    """
    Generate joint distribution statistics for TSO-hamming+G and rcTSO-hamming+C,
    using only TSO/core_TSO for hamming and last 3nt for G/C count.
    Also output grouped percentage curve TSV for G and C mode.
    """
    bam = pysam.AlignmentFile(bam_path, "rb")
    fasta = Fasta(fasta_path)
    seq_len = len(tso_seq)
    core_tso, g3 = get_core_tso(tso_seq)
    corelen = len(core_tso)
    rc_core_tso = reverse_complement(core_tso)
    stat_g = Counter()
    stat_c = Counter()
    logo_records = []
    for read in bam.fetch(until_eof=True):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        chrom = bam.get_reference_name(read.reference_id)
        pos = read.reference_start
        strand = '-' if read.is_reverse else '+'
        seq = fetch_upstream_seq(fasta, chrom, pos, strand, seq_len)
        if len(seq) != seq_len:
            continue
        g_count = count_last3g(seq)
        c_count = count_last3c(seq)
        core_window = seq[:corelen]
        ham_g = hamming_distance_iupac(core_window, core_tso)
        ham_c = hamming_distance_iupac(core_window, rc_core_tso)
        stat_g[(g_count, ham_g)] += 1
        stat_c[(c_count, ham_c)] += 1
        logo_records.append({
            'seq': seq, 'core_window': core_window,
            'G_count': g_count, 'C_count': c_count,
            'hamming_TSO': ham_g, 'hamming_rcTSO': ham_c
        })
    # Output G/C distributions
    pd.DataFrame([(g, h, n) for (g, h), n in stat_g.items()],
                 columns=['G_count', 'Hamming_TSO', 'Count']).to_csv(g_out, sep='\t', index=False)
    pd.DataFrame([(c, h, n) for (c, h), n in stat_c.items()],
                 columns=['C_count', 'Hamming_rcTSO', 'Count']).to_csv(c_out, sep='\t', index=False)
    pd.DataFrame(logo_records).to_csv(logo_table, sep='\t', index=False)
    print(f"Stat outputs: {g_out}, {c_out}, {logo_table}")
    # ---- Output curve tables ----
    total_reads = len(logo_records)
    g_curve = []
    for max_h in range(corelen+1):
        for max_non_g in range(4):
            n = sum(
                (rec['hamming_TSO'] <= max_h) and ((3 - rec['G_count']) <= max_non_g)
                for rec in logo_records
            )
            g_curve.append({'max_hamming': max_h, 'max_non_g': max_non_g, 'percent_reads': 100*n/total_reads})
    g_curve_df = pd.DataFrame(g_curve)
    if g_curve_out:
        g_curve_df.to_csv(g_curve_out, sep='\t', index=False)
        print(f"G curve saved: {g_curve_out}")
    c_curve = []
    for max_h in range(corelen+1):
        for max_non_c in range(4):
            n = sum(
                (rec['hamming_rcTSO'] <= max_h) and ((3 - rec['C_count']) <= max_non_c)
                for rec in logo_records
            )
            c_curve.append({'max_hamming': max_h, 'max_non_c': max_non_c, 'percent_reads': 100*n/total_reads})
    c_curve_df = pd.DataFrame(c_curve)
    if c_curve_out:
        c_curve_df.to_csv(c_curve_out, sep='\t', index=False)
        print(f"C curve saved: {c_curve_out}")

def plot_mode(g_out, c_out, logo_table, logo_mode='G', logo_png=None, heatmap_png=None,
              hamming_max=None, min_g=0, min_c=0, y_mode='rel',
              g_curve=None, c_curve=None, curve_png=None):
    """
    Visualize stats (heatmap, sequence logo, and percentage curve plots).
    curve_png: output path for grouped percentage curve png.
    """
    # Heatmap
    if g_out:
        df_g = pd.read_csv(g_out, sep='\t')
        df_g = df_g[df_g['G_count'].isin([0,1,2,3])]
        pt = df_g.pivot(index='G_count', columns='Hamming_TSO', values='Count').fillna(0)
        plt.figure(figsize=(6,4))
        plt.title("TSO Hamming vs last 3 G count")
        plt.xlabel("Hamming to TSO core")
        plt.ylabel("G in last 3 nt")
        plt.imshow(pt.values, aspect='auto', origin='lower', extent=[pt.columns.min()-0.5, pt.columns.max()+0.5, -0.5, 3.5])
        plt.colorbar(label="Count")
        plt.yticks([0,1,2,3])
        if heatmap_png:
            plt.savefig(heatmap_png.replace('.png','_g.png'))
        else:
            plt.show()
    if c_out:
        df_c = pd.read_csv(c_out, sep='\t')
        df_c = df_c[df_c['C_count'].isin([0,1,2,3])]
        pt = df_c.pivot(index='C_count', columns='Hamming_rcTSO', values='Count').fillna(0)
        plt.figure(figsize=(6,4))
        plt.title("rcTSO Hamming vs last 3 C count")
        plt.xlabel("Hamming to rcTSO core")
        plt.ylabel("C in last 3 nt")
        plt.imshow(pt.values, aspect='auto', origin='lower', extent=[pt.columns.min()-0.5, pt.columns.max()+0.5, -0.5, 3.5])
        plt.colorbar(label="Count")
        plt.yticks([0,1,2,3])
        if heatmap_png:
            plt.savefig(heatmap_png.replace('.png','_c.png'))
        else:
            plt.show()
    # Sequence logo
    if logo_mode and logo_png:
        df = pd.read_csv(logo_table, sep='\t')
        if logo_mode == 'G':
            df['G_count'] = pd.to_numeric(df['G_count'], errors='coerce')
            sel = (df['G_count'] >= min_g)
            if hamming_max is not None:
                sel &= (df['hamming_TSO'] <= hamming_max)
            seqs = df[sel]['seq'].tolist()
            print(f"hamming_max: {hamming_max}, min_g: {min_g}")
        elif logo_mode == 'C':
            sel = (df['C_count'] >= min_c)
            if hamming_max is not None:
                sel &= (df['hamming_rcTSO'] <= hamming_max)
            seqs = df[sel]['seq'].tolist()
            print(f"hamming_max: {hamming_max}, min_g: {min_c}")
        else:
            seqs = df['seq'].tolist()
        if not seqs:
            print("No sequences pass filtering for logo plot.")
            return
        counts_matrix = logomaker.alignment_to_matrix(seqs)
        if y_mode == 'rel':
            counts_matrix = counts_matrix.div(counts_matrix.sum(axis=1), axis=0)
        elif y_mode == 'bits':
            # logomaker handles pseudocounts and info content, but must use counts
            # Fill NA with 0 (rare edge case)
            counts_matrix = counts_matrix.fillna(0)
            info_matrix = logomaker.transform_matrix(counts_matrix, from_type='counts', to_type='information')
            counts_matrix = info_matrix
        # No transform for y_mode == 'count'
        plt.figure(figsize=(min(15, len(seqs[0])//1.2), 3))
        logomaker.Logo(counts_matrix)
        plt.title(f"Sequence logo (mode {logo_mode}, n={len(seqs)}, y={y_mode})")
        plt.tight_layout()
        if y_mode == 'rel':
            plt.ylim(0, 1.05)
            plt.yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0])
            plt.ylabel('Relative Frequency')
        elif y_mode == 'count':
            plt.ylabel('Count')
        elif y_mode == 'bits':
            plt.ylim(0, 2.1)  # Max bits per base for DNA=2
            plt.ylabel('Bits')
        plt.savefig(logo_png)
        print(f"Total reads in logo_table: {len(df)}")
        print(f"Reads used for logo (after filter): {len(seqs)}")
        print(f"Logo saved: {logo_png}")
    # Percentage curve plot
    if g_curve and logo_mode == 'G':
        df = pd.read_csv(g_curve, sep='\t')
        fig, axs = plt.subplots(1, 4, figsize=(16, 4), sharey=True)
        for idx, ng in enumerate([0, 1, 2, 3]):
            sub = df[df['max_non_g'] == ng]
            axs[idx].plot(sub['max_hamming'], sub['percent_reads'], marker='o')
            axs[idx].set_title(f'≤{ng} upstream non-Gs')
            #axs[idx].axvline(x=5, color='k', ls='--')
            axs[idx].set_xlabel('Max hamming distance')
            if idx == 0:
                axs[idx].set_ylabel('Percent reads [%]')
        plt.tight_layout()
        if curve_png:
            plt.savefig(curve_png)
            print(f"Curve plot saved: {curve_png}")
        else:
            plt.show()
    if c_curve and logo_mode == 'C':
        df = pd.read_csv(c_curve, sep='\t')
        fig, axs = plt.subplots(1, 4, figsize=(16, 4), sharey=True)
        for idx, nc in enumerate([0, 1, 2, 3]):
            sub = df[df['max_non_c'] == nc]
            axs[idx].plot(sub['max_hamming'], sub['percent_reads'], marker='o')
            axs[idx].set_title(f'≤{nc} upstream non-Cs')
            #axs[idx].axvline(x=5, color='k', ls='--')
            axs[idx].set_xlabel('Max hamming distance')
            if idx == 0:
                axs[idx].set_ylabel('Percent reads [%]')
        plt.tight_layout()
        if curve_png:
            plt.savefig(curve_png)
            print(f"Curve plot saved: {curve_png}")
        else:
            plt.show()

def filter_mode(bam_path, fasta_path, tso_seq, mode, min_g=None, min_c=None, max_hamming=None, outbam=None):
    """
    Filters input BAM, writing only non-artifact reads to output BAM according to user-specified thresholds.
    Hamming is computed only on TSO core (without last 3nt GGG/CCC).
    """
    bam = pysam.AlignmentFile(bam_path, "rb")
    fasta = Fasta(fasta_path)
    core_tso, _ = get_core_tso(tso_seq)
    rc_core_tso = reverse_complement(core_tso)
    seq_len = len(tso_seq)
    corelen = len(core_tso)
    out_bam = pysam.AlignmentFile(outbam, "wb", template=bam)
    count_written = 0
    count_total = 0
    for read in bam.fetch(until_eof=True):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        chrom = bam.get_reference_name(read.reference_id)
        pos = read.reference_start
        strand = '-' if read.is_reverse else '+'
        seq = fetch_upstream_seq(fasta, chrom, pos, strand, seq_len)
        if len(seq) != seq_len:
            continue
        count_total += 1
        core_window = seq[:corelen]
        if mode == 'strand_invasion':
            g_count = count_last3g(seq)
            ham = hamming_distance_iupac(core_window, core_tso)
            is_artifact = (g_count >= min_g and ham <= max_hamming)
        else:
            c_count = count_last3c(seq)
            ham = hamming_distance_iupac(core_window, rc_core_tso)
            is_artifact = (c_count >= min_c and ham <= max_hamming)
        if not is_artifact:
            out_bam.write(read)
            count_written += 1
    out_bam.close()
    print(f"Filtering complete: {count_written}/{count_total} reads written to {outbam}")

def main():
    """
    Parses arguments and executes corresponding subcommand (stat, plot, strand_invasion, missing_pairing).
    """
    parser = argparse.ArgumentParser(
        description="A pipeline for strand invasion/missing pairing artifact analysis: stat, plot, filter."
    )
    subparsers = parser.add_subparsers(dest='command', required=True)

    # stat
    p_stat = subparsers.add_parser('stat', help="Stat mode: outputs TSV for TSO-hamming+G, rcTSO-hamming+C, logo table, and grouped percentage curves.")
    p_stat.add_argument('-b', '--bam', required=True)
    p_stat.add_argument('-f', '--fasta', required=True)
    p_stat.add_argument('-t', '--tso_seq', required=True)
    p_stat.add_argument('--g_out', required=True, help='TSV for TSO-hamming+G')
    p_stat.add_argument('--c_out', required=True, help='TSV for rcTSO-hamming+C')
    p_stat.add_argument('--logo_table', required=True, help='TSV for all window sequences and features')
    p_stat.add_argument('--g_curve_out', default=None, help='Grouped percent curve table for G mode')
    p_stat.add_argument('--c_curve_out', default=None, help='Grouped percent curve table for C mode')

    # plot
    p_plot = subparsers.add_parser('plot', help="Plot mode: visualize heatmaps/logo/curve from stat .tsv")
    p_plot.add_argument('--g_out', help='TSV for TSO-hamming+G')
    p_plot.add_argument('--c_out', help='TSV for rcTSO-hamming+C')
    p_plot.add_argument('--logo_table', help='TSV for all window sequences and features')
    p_plot.add_argument('--logo_mode', choices=['G', 'C'], default='G', help='Which mode to plot logo')
    p_plot.add_argument('--logo_png', help='Output path for logo PNG')
    p_plot.add_argument('--heatmap_png', help='Prefix for saving heatmaps (will append _g/_c)')
    p_plot.add_argument('--hamming_max', type=int, default=None, help='Filter: max hamming distance for logo')
    p_plot.add_argument('--min_g', type=int, default=0, help='Filter: min G for logo (if G mode)')
    p_plot.add_argument('--min_c', type=int, default=0, help='Filter: min C for logo (if C mode)')
    p_plot.add_argument('--y_mode', choices=['rel', 'count', 'bits'], default='rel',
    help='Logo Y axis: rel=relative freq, count=absolute count, bits=information content')
    p_plot.add_argument('--g_curve', default=None, help='Grouped percent curve table for G mode')
    p_plot.add_argument('--c_curve', default=None, help='Grouped percent curve table for C mode')
    p_plot.add_argument('--curve_png', default=None, help='Output PNG for curve plot')

    # strand_invasion
    p_g = subparsers.add_parser('strand_invasion', help="Filter by G in last 3 nt + TSO hamming")
    p_g.add_argument('-b', '--bam', required=True)
    p_g.add_argument('-f', '--fasta', required=True)
    p_g.add_argument('-t', '--tso_seq', required=True)
    p_g.add_argument('--min_g', type=int, required=True, help='Minimum G in last 3 nt')
    p_g.add_argument('--max_hamming', type=int, required=True, help='Max hamming distance (TSO core)')
    p_g.add_argument('--outbam', required=True, help='Output BAM file')

    # missing_pairing
    p_c = subparsers.add_parser('missing_pairing', help="Filter by C in last 3 nt + rcTSO hamming")
    p_c.add_argument('-b', '--bam', required=True)
    p_c.add_argument('-f', '--fasta', required=True)
    p_c.add_argument('-t', '--tso_seq', required=True)
    p_c.add_argument('--min_c', type=int, required=True, help='Minimum C in last 3 nt')
    p_c.add_argument('--max_hamming', type=int, required=True, help='Max hamming distance (rcTSO core)')
    p_c.add_argument('--outbam', required=True, help='Output BAM file')

    args = parser.parse_args()

    if args.command == 'stat':
        stat_mode(args.bam, args.fasta, args.tso_seq, args.g_out, args.c_out, args.logo_table,
                  g_curve_out=args.g_curve_out, c_curve_out=args.c_curve_out)
    elif args.command == 'plot':
        plot_mode(args.g_out, args.c_out, args.logo_table, logo_mode=args.logo_mode, logo_png=args.logo_png,
                  heatmap_png=args.heatmap_png, hamming_max=args.hamming_max, min_g=args.min_g, min_c=args.min_c, y_mode=args.y_mode,
                  g_curve=args.g_curve, c_curve=args.c_curve, curve_png=args.curve_png)
    elif args.command == 'strand_invasion':
        filter_mode(args.bam, args.fasta, args.tso_seq, 'strand_invasion', min_g=args.min_g, max_hamming=args.max_hamming, outbam=args.outbam)
    elif args.command == 'missing_pairing':
        filter_mode(args.bam, args.fasta, args.tso_seq, 'missing_pairing', min_c=args.min_c, max_hamming=args.max_hamming, outbam=args.outbam)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
