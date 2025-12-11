#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
libend_artifact_filter.py

Author: Hui Zhou (original idea and structure) + Assistant (polishing)
Version: 2.1.0
Date: 2025-10-28

Purpose
-------
Pipeline to score, visualize and filter 5' and 3' end artefacts
in RNA-seq libraries.

Top-level workflow
------------------
For each read, at one end (5' or 3') and for one or more motifs
(TSO or oligo-dT primers), we compute two scores:

  - score1
  - score2

A read is considered an end artefact if there exists any motif such that:

  score1 <= score1_thresh  AND  score2 <= score2_thresh

The exact meaning of score1 / score2 depends on the mode:

Mode "5p" (TSO strand invasion)
-------------------------------
  - motif: a single TSO sequence (provided via --tso-seq)
  - window: upstream of genomic 5' end, length = len(TSO)
  - score1 = Hamming distance (IUPAC-aware) between TSO core
             (TSO without trailing GGG) and upstream window core.
  - score2 = nonG in last 3 nt of upstream window = 3 - G_count_last3

Mode "3p" (oligo-dT missing priming)
------------------------------------
  - motifs: multiple oligo-dT primers (from TSV: name, sequence)
  - window: downstream of genomic 3' end, fixed length from --read-window-len
           (default is max oligo length)
  - score1 = barcode_hamming = Hamming distance in the leading non-N region
  - score2 = dT_nonT_excess = number of non-T bases at positions
             where oligo expects 'T'

The pipeline has three subcommands:

  1. score:  compute per-read per-motif scores, 2D histograms and
             threshold grids; all controlled by a single --prefix.
  2. plot:   draw heatmaps, coverage curves and optional sequence logos.
  3. filter: remove artefact reads based on chosen thresholds.

All TSVs written by `score` include a first comment line describing
the meaning of score1/score2 under the selected mode.

Dependencies
------------
pysam, pyfaidx, pandas, numpy, matplotlib, (optional) logomaker

License
-------
MIT
"""

from __future__ import annotations
import argparse
import logging
import sys
from typing import Dict, List, Tuple, Optional, Iterable, Any
from collections import Counter, defaultdict

import pysam
from pyfaidx import Fasta
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

try:
    import logomaker  # type: ignore
    HAS_LOGOMAKER = True
except ImportError:
    HAS_LOGOMAKER = False


# -----------------------------
# IUPAC utilities
# -----------------------------
IUPAC_CODES: Dict[str, set] = {
    'A': {'A'}, 'C': {'C'}, 'G': {'G'}, 'T': {'T'},
    'R': {'A', 'G'}, 'Y': {'C', 'T'}, 'S': {'G', 'C'}, 'W': {'A', 'T'},
    'K': {'G', 'T'}, 'M': {'A', 'C'}, 'B': {'C', 'G', 'T'}, 'D': {'A', 'G', 'T'},
    'H': {'A', 'C', 'T'}, 'V': {'A', 'C', 'G'}, 'N': {'A', 'C', 'G', 'T'}
}

# Complete reverse complement including IUPAC, preserving case
RC_MAP = str.maketrans(
    "ACGTRYSWKMBDHVNacgtryswkmbdhvn",
    "TGCAYRWSMKVHDBNtgcayrwsmkvhdbn"
)


def reverse_complement(seq: str) -> str:
    return seq.translate(RC_MAP)[::-1]


def iupac_match(base1: str, base2: str) -> bool:
    s1 = IUPAC_CODES.get(base1.upper(), set())
    s2 = IUPAC_CODES.get(base2.upper(), set())
    return bool(s1 & s2)


def hamming_distance_iupac(seq1: str, seq2: str) -> int:
    assert len(seq1) == len(seq2)
    return sum(0 if iupac_match(a, b) else 1 for a, b in zip(seq1, seq2))


# -----------------------------
# FASTA/BAM window extraction
# -----------------------------
def clamp_window(chrom_len: int, start: int, end: int) -> Optional[Tuple[int, int]]:
    if start < 0:
        start = 0
    if end > chrom_len:
        end = chrom_len
    if end <= start:
        return None
    return start, end


def fetch_upstream_from_5prime(
    fasta: Fasta,
    chrom: str,
    fivep_pos: int,
    strand: str,
    length: int
) -> Optional[str]:
    """
    Get upstream sequence of given length from a read's 5' mapping coordinate.
    fivep_pos: genomic coordinate (0-based, inclusive) of the 5' end of the read.

    For '+' strand reads:
        upstream is [fivep - L, fivep) on genome (forward strand).
    For '-' strand reads:
        fivep_pos is genomic coordinate of rightmost base;
        upstream in read orientation corresponds to (fivep, fivep + L] on genome,
        which is then reverse-complemented.

    Returns sequence in read 5'->3' orientation (upper-case),
    or None if the window cannot be retrieved with full length.
    """
    chrom_len = len(fasta[chrom])
    if strand == '+':
        w = clamp_window(chrom_len, fivep_pos - length, fivep_pos)
        if w is None or (w[1] - w[0]) != length:
            return None
        return str(fasta[chrom][w[0]:w[1]].seq).upper()
    else:
        w = clamp_window(chrom_len, fivep_pos + 1, fivep_pos + 1 + length)
        if w is None or (w[1] - w[0]) != length:
            return None
        seq = str(fasta[chrom][w[0]:w[1]].seq)
        return reverse_complement(seq).upper()


def fetch_downstream_from_3prime(
    fasta: Fasta,
    chrom: str,
    threep_pos: int,
    strand: str,
    length: int
) -> Optional[str]:
    """
    Get downstream sequence of given length from a read's 3' mapping coordinate.
    threep_pos: genomic coordinate (0-based, inclusive) of the 3' end of the read.

    For '+' strand reads:
        downstream is (threep, threep + L] on genome.
    For '-' strand reads:
        downstream in read orientation corresponds to [threep - L + 1, threep] on
        genome (leftwards), which is then reverse-complemented.

    Returns sequence in read 5'->3' orientation (upper-case),
    or None if the window cannot be retrieved with full length.
    """
    chrom_len = len(fasta[chrom])
    if strand == '+':
        w = clamp_window(chrom_len, threep_pos + 1, threep_pos + 1 + length)
        if w is None or (w[1] - w[0]) != length:
            return None
        return str(fasta[chrom][w[0]:w[1]].seq).upper()
    else:
        w = clamp_window(chrom_len, threep_pos - length + 1, threep_pos + 1)
        if w is None or (w[1] - w[0]) != length:
            return None
        seq = str(fasta[chrom][w[0]:w[1]].seq)
        return reverse_complement(seq).upper()


# -----------------------------
# 5' mode helpers (TSO)
# -----------------------------
def get_core_tso(tsoseq: str, g3_len: int = 3) -> Tuple[str, str]:
    tso = tsoseq.upper()
    assert len(tso) > g3_len, "TSO length must be > g3_len"
    core = tso[:-g3_len]
    tail = tso[-g3_len:]
    return core, tail


def count_last3(seq: str, base: str) -> int:
    return seq[-3:].count(base)


# -----------------------------
# 3' mode helpers (oligos)
# -----------------------------
class OligoSpec:
    """
    Represents an oligo-dT primer with masks:
      - barcode_mask: indices to evaluate Hamming against oligo's bases
      - t_mask: indices where oligo expects 'T'
      - total_len: length of oligo sequence template
    """

    def __init__(self, name: str, sequence: str):
        self.name = name
        self.seq = sequence.upper()
        self.total_len = len(self.seq)
        self.barcode_mask = self._infer_barcode_mask()
        self.t_mask = self._infer_t_mask()
        self._validate()

    def _infer_barcode_mask(self) -> List[int]:
        # Barcode: leading consecutive non-N characters from 5' end
        mask: List[int] = []
        for i, ch in enumerate(self.seq):
            if ch == 'N':
                break
            mask.append(i)
        return mask

    def _infer_t_mask(self) -> List[int]:
        return [i for i, ch in enumerate(self.seq) if ch == 'T']

    def _validate(self) -> None:
        if self.total_len <= 0:
            raise ValueError(f"Oligo {self.name}: empty sequence.")
        if any(i < 0 or i >= self.total_len for i in self.barcode_mask):
            raise ValueError(f"Oligo {self.name}: barcode indices out of range.")
        if any(i < 0 or i >= self.total_len for i in self.t_mask):
            raise ValueError(f"Oligo {self.name}: T-mask indices out of range.")
        if len(self.barcode_mask) == 0:
            logging.warning(
                "Oligo %s: no leading barcode (no non-N before first N).",
                self.name
            )

    def barcode_hamming(self, window_seq: str) -> int:
        if not self.barcode_mask:
            return 0
        if len(window_seq) < self.total_len:
            return np.inf
        mismatches = 0
        for i in self.barcode_mask:
            if window_seq[i] != self.seq[i]:
                mismatches += 1
        return mismatches

    def dt_nonT_excess(self, window_seq: str) -> int:
        if len(window_seq) < self.total_len:
            return np.inf
        cnt = 0
        for i in self.t_mask:
            if window_seq[i] != 'T':
                cnt += 1
        return cnt


def load_oligos(tsv_path: str) -> List[OligoSpec]:
    import re
    df = pd.read_csv(tsv_path, sep="\t")
    lower_cols = {c.lower(): c for c in df.columns}
    if not {'name', 'sequence'}.issubset(set(lower_cols.keys())):
        raise ValueError("Oligo TSV must have columns: name, sequence (tab-delimited).")
    name_col = lower_cols['name']
    seq_col = lower_cols['sequence']

    valid = set("ACGTRYSWKMBDHVN")
    oligos: List[OligoSpec] = []
    for _, row in df.iterrows():
        name = str(row[name_col]).strip()
        raw = str(row[seq_col])
        seq = re.sub(r"\s+", "", raw).upper()
        seq = "".join(ch for ch in seq if ch in valid)
        if not seq:
            raise ValueError(f"Oligo {name}: empty after cleaning. Check your file.")
        oligos.append(OligoSpec(name, seq))
        logging.info("Loaded oligo %s: len=%d (cleaned).", name, len(seq))
    return oligos


# -----------------------------
# Generic read iterator
# -----------------------------
def iterate_reads(
    bam: pysam.AlignmentFile,
    min_mapq: int = 0,
    max_reads: Optional[int] = None
) -> Iterable[pysam.AlignedSegment]:
    n = 0
    for read in bam.fetch(until_eof=True):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        if read.mapping_quality < min_mapq:
            continue
        yield read
        n += 1
        if max_reads is not None and n >= max_reads:
            break


# -----------------------------
# Score semantics helper
# -----------------------------
def get_score_definitions(mode: str) -> Tuple[str, str]:
    if mode == "5p":
        return (
            "score1 = Hamming distance to TSO core (IUPAC-aware)",
            "score2 = non-G count in last 3 nt (3 - G_count_last3)"
        )
    else:
        return (
            "score1 = barcode_hamming (Hamming in leading non-N barcode region)",
            "score2 = dT_nonT_excess (non-T count at oligo 'T' positions)"
        )


def get_score_axis_labels(mode: str) -> Tuple[str, str]:
    if mode == "5p":
        return (
            "score1 (Hamming to TSO core)",
            "score2 (non-G in last 3 nt)"
        )
    else:
        return (
            "score1 (barcode Hamming)",
            "score2 (dT non-T excess)"
        )


# -----------------------------
# Per-read scoring functions
# -----------------------------
def score_read_5p(
    read: pysam.AlignedSegment,
    bam: pysam.AlignmentFile,
    fasta: Fasta,
    tso_seq: str
) -> Optional[Dict[str, Any]]:
    chrom = bam.get_reference_name(read.reference_id)
    strand = '-' if read.is_reverse else '+'
    fivep = (read.reference_end - 1) if read.is_reverse else read.reference_start
    if chrom not in fasta:
        return None

    seq_len = len(tso_seq)
    core_tso, _ = get_core_tso(tso_seq, g3_len=3)
    corelen = len(core_tso)

    window = fetch_upstream_from_5prime(fasta, chrom, fivep, strand, seq_len)
    if window is None:
        return None

    g_count = count_last3(window, 'G')
    non_g = 3 - g_count
    core_window = window[:corelen]
    ham = hamming_distance_iupac(core_window, core_tso)

    return {
        "read_qname": read.query_name,
        "chrom": chrom,
        "strand": strand,
        "end_pos": fivep,
        "motif_id": "TSO",
        "score1": int(ham),
        "score2": int(non_g),
        "seq_window": window,
        "g_count_last3": int(g_count)
    }


def score_read_3p(
    read: pysam.AlignedSegment,
    bam: pysam.AlignmentFile,
    fasta: Fasta,
    oligos: List[OligoSpec],
    window_len: int
) -> List[Dict[str, Any]]:
    chrom = bam.get_reference_name(read.reference_id)
    strand = '-' if read.is_reverse else '+'
    threep = (read.reference_start) if read.is_reverse else (read.reference_end - 1)
    if chrom not in fasta:
        return []

    win = fetch_downstream_from_3prime(fasta, chrom, threep, strand, window_len)
    if win is None:
        return []

    rows: List[Dict[str, Any]] = []
    for ol in oligos:
        if len(win) < ol.total_len:
            continue
        w = win[:ol.total_len]
        bh = ol.barcode_hamming(w)
        dt = ol.dt_nonT_excess(w)
        rows.append({
            "read_qname": read.query_name,
            "chrom": chrom,
            "strand": strand,
            "end_pos": threep,
            "motif_id": ol.name,
            "score1": int(bh),
            "score2": int(dt),
            "seq_window": w,
            "barcode_hamming_raw": int(bh),
            "dT_nonT_excess_raw": int(dt)
        })
    return rows


# -----------------------------
# SCORE subcommand
# -----------------------------
def score_mode(
    mode: str,
    bam_path: str,
    fasta_path: str,
    tso_seq: Optional[str],
    oligo_tsv: Optional[str],
    read_window_len: Optional[int],
    prefix: str,
    score1_max_grid: int,
    score2_max_grid: int,
    min_mapq: int,
    max_reads: Optional[int]
) -> None:
    logging.info("Score: starting (mode=%s)", mode)
    bam = pysam.AlignmentFile(bam_path, "rb")
    fasta = Fasta(fasta_path)

    oligos: Optional[List[OligoSpec]] = None
    if mode == "3p":
        if oligo_tsv is None:
            raise ValueError("Mode 3p requires --oligos TSV.")
        oligos = load_oligos(oligo_tsv)
        if not oligos:
            raise ValueError("No oligos loaded.")
        if read_window_len is None:
            read_window_len = max(o.total_len for o in oligos)
        logging.info("3p: using read_window_len=%d", read_window_len)
    else:
        if tso_seq is None:
            raise ValueError("Mode 5p requires --tso-seq.")
        tso_seq = tso_seq.upper()

    per_rows: List[Dict[str, Any]] = []
    hist_counter: Counter[Tuple[int, int]] = Counter()
    pairs_by_read: defaultdict[str, List[Tuple[int, int]]] = defaultdict(list)

    total_reads = 0
    scored_reads = 0
    skipped_window = 0

    for read in iterate_reads(bam, min_mapq=min_mapq, max_reads=max_reads):
        total_reads += 1
        if mode == "5p":
            assert tso_seq is not None
            row = score_read_5p(read, bam, fasta, tso_seq)
            if row is None:
                skipped_window += 1
                continue
            per_rows.append(row)
            scored_reads += 1
            s1 = int(row["score1"])
            s2 = int(row["score2"])
            hist_counter[(s1, s2)] += 1
            pairs_by_read[row["read_qname"]].append((s1, s2))
        else:
            assert oligos is not None
            assert read_window_len is not None
            rows = score_read_3p(read, bam, fasta, oligos, read_window_len)
            if not rows:
                skipped_window += 1
                continue
            scored_reads += 1
            for r in rows:
                per_rows.append(r)
                s1 = int(r["score1"])
                s2 = int(r["score2"])
                hist_counter[(s1, s2)] += 1
                pairs_by_read[r["read_qname"]].append((s1, s2))

    bam.close()
    fasta.close()

    logging.info(
        "Score: total_reads=%d, scored_reads=%d, skipped_window=%d",
        total_reads, scored_reads, skipped_window
    )

    if not per_rows:
        logging.warning("No scores computed; no output will be written.")
        return

    score1_desc, score2_desc = get_score_definitions(mode)

    # 1) per-read table
    per_df = pd.DataFrame(per_rows)

    per_path = f"{prefix}.score_per_read.tsv"
    with open(per_path, "w") as fh:
        fh.write(
            f"# mode={mode}, {score1_desc}; {score2_desc}\n"
        )
        per_df.to_csv(fh, sep="\t", index=False)
    logging.info("Score: per-read table written to %s", per_path)

    # 2) 2D histogram
    hist_records = [(s1, s2, n) for (s1, s2), n in hist_counter.items()]
    hist_df = pd.DataFrame(hist_records, columns=["score1", "score2", "count"])
    hist_path = f"{prefix}.score_2d_hist.tsv"
    with open(hist_path, "w") as fh:
        fh.write(
            f"# mode={mode}, {score1_desc}; {score2_desc}\n"
        )
        hist_df.to_csv(fh, sep="\t", index=False)
    logging.info("Score: 2D histogram written to %s", hist_path)

    # 3) threshold grid (ANY-motif condition per read)
    all_s1 = [s1 for (s1, _s2) in hist_counter.keys()]
    all_s2 = [s2 for (_s1, s2) in hist_counter.keys()]
    max_s1_obs = max(all_s1)
    max_s2_obs = max(all_s2)
    s1_cap = min(score1_max_grid, max_s1_obs)
    s2_cap = min(score2_max_grid, max_s2_obs)
    logging.info(
        "Score: grid caps score1_max=%d, score2_max=%d (observed max_s1=%d, max_s2=%d)",
        s1_cap, s2_cap, max_s1_obs, max_s2_obs
    )

    unique_reads = list(pairs_by_read.keys())
    n_reads = len(unique_reads)
    grid_records: List[Dict[str, Any]] = []

    for s1_max in range(s1_cap + 1):
        for s2_max in range(s2_cap + 1):
            flagged = 0
            for rq in unique_reads:
                plist = pairs_by_read.get(rq, [])
                if any((s1 <= s1_max and s2 <= s2_max) for (s1, s2) in plist):
                    flagged += 1
            pct = 100.0 * flagged / n_reads if n_reads else 0.0
            grid_records.append({
                "score1_max": s1_max,
                "score2_max": s2_max,
                "percent_reads_flagged": pct
            })

    grid_df = pd.DataFrame(grid_records)
    grid_path = f"{prefix}.threshold_grid.tsv"
    with open(grid_path, "w") as fh:
        fh.write(
            f"# mode={mode}, {score1_desc}; {score2_desc}\n"
        )
        grid_df.to_csv(fh, sep="\t", index=False)
    logging.info("Score: threshold grid written to %s", grid_path)


# -----------------------------
# PLOT subcommand
# -----------------------------
def plot_mode(
    mode: str,
    prefix: str,
    logo_mode: Optional[str],
    logo_score1_max: Optional[int],
    logo_score2_max: Optional[int],
    group_by_motif: bool
) -> None:
    """
    logo_mode:
        None  -> no logo
        'rel' -> relative frequency
        'bits'/'bit' -> information content (bits)
        'count' -> raw counts
    """
    logging.info("Plot: starting (mode=%s)", mode)

    hist_path = f"{prefix}.score_2d_hist.tsv"
    grid_path = f"{prefix}.threshold_grid.tsv"
    per_path = f"{prefix}.score_per_read.tsv"

    def _read_tsv(path: str) -> pd.DataFrame:
        return pd.read_csv(path, sep="\t", comment="#")

    # 1) 2D histogram heatmap
    try:
        hist_df = _read_tsv(hist_path)
    except FileNotFoundError:
        logging.error("Histogram TSV not found: %s", hist_path)
        hist_df = None

    x_label, y_label = get_score_axis_labels(mode)
    score1_desc, score2_desc = get_score_definitions(mode)

    if hist_df is not None and not hist_df.empty:
        pt = hist_df.pivot(index="score2", columns="score1", values="count").fillna(0)
        plt.figure(figsize=(7, 5))
        im = plt.imshow(
            pt.values,
            aspect="auto",
            origin="lower",
            extent=[
                pt.columns.min() - 0.5,
                pt.columns.max() + 0.5,
                pt.index.min() - 0.5,
                pt.index.max() + 0.5
            ]
        )
        plt.colorbar(im, label="Count")
        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(f"{mode} mode: 2D histogram of scores\n{score1_desc}; {score2_desc}")
        plt.tight_layout()
        out_png = f"{prefix}.score2d_heatmap.png"
        plt.savefig(out_png, dpi=200, bbox_inches="tight")
        plt.close()
        logging.info("Plot: 2D histogram heatmap saved to %s", out_png)
    else:
        logging.warning("No histogram data to plot.")

    # 2) threshold grid heatmap + curves
    try:
        grid_df = _read_tsv(grid_path)
    except FileNotFoundError:
        logging.error("Grid TSV not found: %s", grid_path)
        grid_df = None

    if grid_df is not None and not grid_df.empty:
        ptg = grid_df.pivot(
            index="score2_max",
            columns="score1_max",
            values="percent_reads_flagged"
        ).fillna(0)
        plt.figure(figsize=(7, 5))
        im = plt.imshow(
            ptg.values,
            aspect="auto",
            origin="lower",
            extent=[
                ptg.columns.min() - 0.5,
                ptg.columns.max() + 0.5,
                ptg.index.min() - 0.5,
                ptg.index.max() + 0.5
            ],
            vmin=0,
            vmax=100.0
        )
        plt.colorbar(im, label="Percent reads flagged [%]")
        plt.xlabel("score1_max")
        plt.ylabel("score2_max")
        plt.title(
            f"{mode} mode: coverage over thresholds (ANY motif)\n"
            f"{score1_desc}; {score2_desc}"
        )
        from matplotlib.ticker import MaxNLocator
        ax = plt.gca()
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        plt.tight_layout()
        out_png = f"{prefix}.threshold_grid_heatmap.png"
        plt.savefig(out_png, dpi=200, bbox_inches="tight")
        plt.close()
        logging.info("Plot: threshold grid heatmap saved to %s", out_png)

        # coverage curves
        unique_s2 = sorted(grid_df["score2_max"].unique())
        if len(unique_s2) <= 6:
            slice_s2 = unique_s2
        else:
            idxs = [
                0,
                len(unique_s2) // 5,
                (2 * len(unique_s2)) // 5,
                (3 * len(unique_s2)) // 5,
                (4 * len(unique_s2)) // 5,
                len(unique_s2) - 1
            ]
            slice_s2 = [unique_s2[i] for i in sorted(set(idxs))]

        plt.figure(figsize=(8, 5))
        for s2 in slice_s2:
            sub = grid_df[grid_df["score2_max"] == s2].sort_values("score1_max")
            plt.plot(
                sub["score1_max"],
                sub["percent_reads_flagged"],
                marker="o",
                label=f"score2≤{int(s2)}"
            )
        plt.xlabel("score1_max")
        plt.ylabel("Percent reads flagged [%]")
        plt.title(
            f"{mode} mode: coverage curves\n"
            f"{score1_desc}; {score2_desc}"
        )
        from matplotlib.ticker import MaxNLocator as MNL
        ax = plt.gca()
        ax.xaxis.set_major_locator(MNL(integer=True))
        plt.legend()
        plt.tight_layout()
        out_png = f"{prefix}.threshold_slices.png"
        plt.savefig(out_png, dpi=200, bbox_inches="tight")
        plt.close()
        logging.info("Plot: threshold coverage curves saved to %s", out_png)
    else:
        logging.warning("No threshold grid data to plot.")

    # 3) sequence logo
    if logo_mode is None:
        return

    if not HAS_LOGOMAKER:
        logging.error("logomaker not installed; logo plot skipped.")
        return

    # 支持 'bit' 作为 'bits' 的别名
    if logo_mode == "bit":
        logo_mode = "bits"

    if logo_mode not in {"rel", "bits", "count"}:
        logging.error("Unknown logo mode: %s (use rel/bits/count)", logo_mode)
        return

    if logo_score1_max is None or logo_score2_max is None:
        logging.error(
            "Logo requested but --logo-score1-max and --logo-score2-max are required."
        )
        return

    try:
        per_df = _read_tsv(per_path)
    except FileNotFoundError:
        logging.error("Per-read TSV not found: %s", per_path)
        return

    if "seq_window" not in per_df.columns:
        logging.error("Per-read TSV has no 'seq_window' column; cannot plot logo.")
        return

    mask = (
        (per_df["score1"] <= logo_score1_max) &
        (per_df["score2"] <= logo_score2_max)
    )
    df_sel = per_df.loc[mask].copy()
    if df_sel.empty:
        logging.warning(
            "No reads meet logo criteria (score1≤%d, score2≤%d).",
            logo_score1_max, logo_score2_max
        )
        return

    df_sel["seq_window"] = df_sel["seq_window"].astype(str)
    lengths = df_sel["seq_window"].map(len)
    length_mode = int(lengths.mode().iloc[0])
    if not (lengths == length_mode).all():
        bad = (lengths != length_mode).sum()
        logging.info(
            "Logo: %d rows have seq_window length != %d and will be ignored.",
            bad, length_mode
        )
        df_sel = df_sel[lengths == length_mode]
        if df_sel.empty:
            logging.warning("No sequences left for logo after length harmonization.")
            return

    def _draw_logo(seqs: List[str], out_png: str, title_suffix: str) -> None:
        cm = logomaker.alignment_to_matrix(seqs)
        y_label = ""
        if logo_mode == "rel":
            cm = cm.div(cm.sum(axis=1), axis=0)
            y_label = "Relative frequency"
        elif logo_mode == "bits":
            cm = cm.fillna(0)
            cm = logomaker.transform_matrix(cm, from_type="counts", to_type="information")
            y_label = "Bits"
        else:  # 'count'
            y_label = "Count"

        plt.figure(figsize=(min(15, len(seqs[0]) / 1.2), 3))
        logomaker.Logo(cm)
        plt.title(
            f"{mode} mode logo {title_suffix}\n"
            f"{score1_desc}; {score2_desc}\n"
            f"(n={len(seqs)}, score1≤{logo_score1_max}, score2≤{logo_score2_max}, y={logo_mode})"
        )
        if logo_mode == "rel":
            plt.ylim(0, 1.05)
        elif logo_mode == "bits":
            plt.ylim(0, 2.1)
        plt.ylabel(y_label)
        plt.tight_layout()
        plt.savefig(out_png, dpi=200, bbox_inches="tight")
        plt.close()
        logging.info("Logo saved: %s (n=%d, mode=%s)", out_png, len(seqs), logo_mode)

    # all combined
    seqs_all = df_sel["seq_window"].astype(str).tolist()
    _draw_logo(seqs_all, f"{prefix}.logo_all.png", "(all motifs)")

    # optionally group by motif
    if group_by_motif and "motif_id" in df_sel.columns:
        base = f"{prefix}.logo"
        ext = ".png"
        for motif, sub in df_sel.groupby("motif_id", dropna=False):
            seqs = sub["seq_window"].astype(str).tolist()
            safe = str(motif).replace("/", "_").replace("\\", "_").replace(" ", "_")
            out_png = f"{base}_{safe}{ext}"
            _draw_logo(seqs, out_png, f"(motif={safe})")


# -----------------------------
# FILTER subcommand
# -----------------------------
def filter_mode(
    mode: str,
    bam_path: str,
    fasta_path: str,
    tso_seq: Optional[str],
    oligo_tsv: Optional[str],
    read_window_len: Optional[int],
    outbam: str,
    score1_thresh: int,
    score2_thresh: int,
    min_mapq: int,
    max_reads: Optional[int]
) -> None:
    logging.info(
        "Filter: starting (mode=%s, score1_thresh=%d, score2_thresh=%d)",
        mode, score1_thresh, score2_thresh
    )
    bam = pysam.AlignmentFile(bam_path, "rb")
    out_bam = pysam.AlignmentFile(outbam, "wb", template=bam)
    fasta = Fasta(fasta_path)

    oligos: Optional[List[OligoSpec]] = None
    if mode == "3p":
        if oligo_tsv is None:
            raise ValueError("Mode 3p requires --oligos TSV.")
        oligos = load_oligos(oligo_tsv)
        if not oligos:
            raise ValueError("No oligos loaded.")
        if read_window_len is None:
            read_window_len = max(o.total_len for o in oligos)
        logging.info("3p: using read_window_len=%d", read_window_len)
    else:
        if tso_seq is None:
            raise ValueError("Mode 5p requires --tso-seq.")
        tso_seq = tso_seq.upper()

    total = 0
    written = 0
    artefact = 0
    skipped_window = 0

    for read in iterate_reads(bam, min_mapq=min_mapq, max_reads=max_reads):
        total += 1
        scores: List[Tuple[int, int]] = []

        if mode == "5p":
            assert tso_seq is not None
            row = score_read_5p(read, bam, fasta, tso_seq)
            if row is None:
                skipped_window += 1
            else:
                scores.append((int(row["score1"]), int(row["score2"])))
        else:
            assert oligos is not None
            assert read_window_len is not None
            rows = score_read_3p(read, bam, fasta, oligos, read_window_len)
            if not rows:
                skipped_window += 1
            else:
                for r in rows:
                    scores.append((int(r["score1"]), int(r["score2"])))

        # conservative: if we could not score this read, keep it
        if not scores:
            out_bam.write(read)
            written += 1
            continue

        is_art = any(
            (s1 <= score1_thresh and s2 <= score2_thresh)
            for (s1, s2) in scores
        )
        if is_art:
            artefact += 1
        else:
            out_bam.write(read)
            written += 1

    bam.close()
    out_bam.close()
    fasta.close()

    logging.info(
        "Filter: total_reads=%d, written=%d, artefact=%d, skipped_window=%d",
        total, written, artefact, skipped_window
    )


# -----------------------------
# CLI
# -----------------------------
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Unified 5' (TSO) & 3' (oligo-dT) artefact analysis: score, plot, filter."
    )
    p.add_argument(
        "--mode",
        choices=["5p", "3p"],
        required=True,
        help="Library end mode: 5p (TSO) or 3p (oligo-dT)"
    )

    sub = p.add_subparsers(dest="cmd", required=True)

    # score
    sp_score = sub.add_parser("score", help="Compute per-read scores and threshold grids.")
    sp_score.add_argument("-b", "--bam", required=True, help="Input BAM")
    sp_score.add_argument("-f", "--fasta", required=True, help="Reference FASTA")
    sp_score.add_argument(
        "--tso-seq",
        help="TSO sequence (ending with GGG). Required for --mode 5p."
    )
    sp_score.add_argument(
        "--oligos",
        help="TSV with oligo primers (name, sequence). Required for --mode 3p."
    )
    sp_score.add_argument(
        "--read-window-len",
        type=int,
        default=None,
        help="3p: window length from 3' end; default = max oligo length."
    )
    sp_score.add_argument(
        "--prefix",
        required=True,
        help="Output prefix for score TSVs."
    )
    sp_score.add_argument(
        "--score1-max-grid",
        type=int,
        default=20,
        help="Max score1 for threshold grid (cap)."
    )
    sp_score.add_argument(
        "--score2-max-grid",
        type=int,
        default=30,
        help="Max score2 for threshold grid (cap)."
    )
    sp_score.add_argument(
        "--min-mapq",
        type=int,
        default=0,
        help="Minimum MAPQ to keep reads."
    )
    sp_score.add_argument(
        "--max-reads",
        type=int,
        default=None,
        help="Optional cap on number of reads processed (for testing)."
    )

    # plot
    sp_plot = sub.add_parser("plot", help="Plot heatmaps, coverage curves and logos.")
    sp_plot.add_argument(
        "--prefix",
        required=True,
        help="Same prefix used in the score step."
    )
    sp_plot.add_argument(
        "--make-logo",
        choices=["rel", "bits", "bit", "count"],
        default=None,
        help="Draw sequence logo with given y-axis mode: rel / bits / count."
    )
    sp_plot.add_argument(
        "--logo-score1-max",
        type=int,
        default=None,
        help="Max score1 for selecting reads used in logo."
    )
    sp_plot.add_argument(
        "--logo-score2-max",
        type=int,
        default=None,
        help="Max score2 for selecting reads used in logo."
    )
    sp_plot.add_argument(
        "--group-by-motif",
        action="store_true",
        help="For logo: draw one logo per motif (oligo/TSO)."
    )

    # filter
    sp_filter = sub.add_parser("filter", help="Filter BAM by artefact rules.")
    sp_filter.add_argument("-b", "--bam", required=True, help="Input BAM")
    sp_filter.add_argument("-f", "--fasta", required=True, help="Reference FASTA")
    sp_filter.add_argument("-o", "--outbam", required=True, help="Output BAM")
    sp_filter.add_argument(
        "--tso-seq",
        help="TSO sequence (ending with GGG). Required for --mode 5p."
    )
    sp_filter.add_argument(
        "--oligos",
        help="TSV with oligo primers (name, sequence). Required for --mode 3p."
    )
    sp_filter.add_argument(
        "--read-window-len",
        type=int,
        default=None,
        help="3p: downstream window length from 3' end."
    )
    sp_filter.add_argument(
        "--score1-thresh",
        type=int,
        required=True,
        help="Threshold for score1 (see mode-specific definition)."
    )
    sp_filter.add_argument(
        "--score2-thresh",
        type=int,
        required=True,
        help="Threshold for score2 (see mode-specific definition)."
    )
    sp_filter.add_argument(
        "--min-mapq",
        type=int,
        default=0,
        help="Minimum MAPQ to keep reads."
    )
    sp_filter.add_argument(
        "--max-reads",
        type=int,
        default=None,
        help="Optional cap on number of reads processed (for testing)."
    )

    # logging
    p.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity."
    )
    return p


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="[%(levelname)s] %(message)s"
    )

    if args.cmd == "score":
        score_mode(
            mode=args.mode,
            bam_path=args.bam,
            fasta_path=args.fasta,
            tso_seq=args.tso_seq,
            oligo_tsv=args.oligos,
            read_window_len=args.read_window_len,
            prefix=args.prefix,
            score1_max_grid=args.score1_max_grid,
            score2_max_grid=args.score2_max_grid,
            min_mapq=args.min_mapq,
            max_reads=args.max_reads
        )

    elif args.cmd == "plot":
        plot_mode(
            mode=args.mode,
            prefix=args.prefix,
            logo_mode=args.make_logo,
            logo_score1_max=args.logo_score1_max,
            logo_score2_max=args.logo_score2_max,
            group_by_motif=args.group_by_motif
        )

    elif args.cmd == "filter":
        filter_mode(
            mode=args.mode,
            bam_path=args.bam,
            fasta_path=args.fasta,
            tso_seq=args.tso_seq,
            oligo_tsv=args.oligos,
            read_window_len=args.read_window_len,
            outbam=args.outbam,
            score1_thresh=args.score1_thresh,
            score2_thresh=args.score2_thresh,
            min_mapq=args.min_mapq,
            max_reads=args.max_reads
        )
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
