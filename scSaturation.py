#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# created by: Hui ZHou
# date: 2025-Jun-19

"""
scSaturation.py

Single-cell sequencing saturation curves for general bam file with libray and per cell support.

Modes (mutually exclusive via --mode):
  - global   : Monte Carlo subsampling of the whole BAM and counting uniques
               --dedup cell       -> unique (CB, UMI)
               --dedup library    -> unique UMI (ignore CB & Gene)
               --dedup cell-gene  -> unique (CB, Gene, UMI)
  - per-cell : For each selected CB, compute EXPECTED unique molecules vs reads
               analytically (no Monte Carlo)
               --per-cell-dedup {cell, cell-gene}

Gene tags:
  --gene-tags "GX,GN" (priority order; default)  |  --keep-missing-gene to keep as "__NA__"

Axes:
  - X axis is always "Number of reads (M)"; you can override both axes via --xlabel/--ylabel
  - Journal style (no grid). Optional --logx, and --no-band (global only)

CSV (long-format):
  - global  : reads,reads_M,mean_unique,std_unique,dedup
  - per-cell: cell,reads,reads_M,expected_unique,dedup

Examples
pdc: ml bioinfo-tools; ml python/3.12.3
--------
# Global, per-cell dedup, 25 points, CSV
python scSaturation.py --mode global --bam in.bam --out-prefix sample \
  --dedup cell --num-points 25 --skip-secondary --csv

# Global, library-level unique UMI (ignore CB/Gene), log-x, no band
python scSaturation.py --mode global --bam in.bam --out-prefix sample_lib \
  --dedup library --logx --no-band

# Global, cell-gene dedup using GX then GN
python scSaturation.py --mode global --bam in.bam --out-prefix sample_cg \
  --dedup cell-gene --gene-tags GX,GN

# Per-cell curves, per-cell dedup, top100 cells (>=2k reads), CSV
python scSaturation.py --mode per-cell --bam in.bam --out-prefix sample_cells \
  --per-cell-dedup cell --per-cell-max 100 --per-cell-min-reads 2000 --csv

# Per-cell curves, cell-gene dedup for a custom CB list
python scSaturation.py --mode per-cell --bam in.bam --out-prefix sample_cells_cg \
  --per-cell-dedup cell-gene --cb-list cells.txt
"""

import argparse
import math
import random
from collections import Counter, defaultdict
from typing import List, Tuple, Dict, Optional

import numpy as np
import pysam
import matplotlib as mpl
import matplotlib.pyplot as plt

try:
    from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    HAS_TQDM = False


# ----------------------------- helpers: tags & I/O ----------------------------

def parse_gene_tags(gene_tags_str: str) -> List[str]:
    toks = gene_tags_str.replace(",", " ").split()
    return [t.strip() for t in toks if t.strip()]


def get_gene_from_read(read: pysam.AlignedSegment, gene_tags: List[str]) -> Optional[str]:
    for tag in gene_tags:
        if read.has_tag(tag):
            g = read.get_tag(tag)
            if g:
                return str(g)
    return None


def load_records(bam_path: str,
                 cb_tag: str,
                 umi_tag: str,
                 gene_tags: List[str],
                 min_mapq: int,
                 skip_secondary: bool,
                 requires_gene: bool,
                 keep_missing_gene: bool) -> List[Tuple[str, str, Optional[str]]]:
    """
    Load records as (CB, UMI, Gene or None). If requires_gene=True:
      - drop reads missing gene unless keep_missing_gene=True (then Gene="__NA__").
    """
    recs: List[Tuple[str, str, Optional[str]]] = []
    with pysam.AlignmentFile(bam_path, "rb") as bam:
        it = bam.fetch(until_eof=True)
        for read in it:
            if skip_secondary and (read.is_secondary or read.is_supplementary):
                continue
            if read.mapping_quality < min_mapq:
                continue
            if not (read.has_tag(cb_tag) and read.has_tag(umi_tag)):
                continue
            cb = read.get_tag(cb_tag)
            ub = read.get_tag(umi_tag)
            if not (cb and ub):
                continue
            gene = get_gene_from_read(read, gene_tags)
            if requires_gene and gene is None:
                if keep_missing_gene:
                    gene = "__NA__"
                else:
                    continue
            recs.append((str(cb), str(ub), None if gene is None else str(gene)))
    return recs


def load_cb_item_counts(bam_path: str,
                        cb_tag: str,
                        umi_tag: str,
                        gene_tags: List[str],
                        min_mapq: int,
                        skip_secondary: bool,
                        per_cell_dedup: str,
                        keep_missing_gene: bool) -> Tuple[Dict[str, Counter], Dict[str, int]]:
    """
    For per-cell curves: return
      cb2items: dict[CB] -> Counter({item_key: multiplicity_in_reads})
                 item_key = UMI (per_cell_dedup='cell')
                          or (Gene, UMI) (per_cell_dedup='cell-gene')
      cb2N    : dict[CB] -> total reads contributing to items for that CB
    """
    cb2items: Dict[str, Counter] = defaultdict(Counter)
    cb2N: Dict[str, int] = defaultdict(int)
    require_gene = (per_cell_dedup == "cell-gene")

    with pysam.AlignmentFile(bam_path, "rb") as bam:
        it = bam.fetch(until_eof=True)
        for read in it:
            if skip_secondary and (read.is_secondary or read.is_supplementary):
                continue
            if read.mapping_quality < min_mapq:
                continue
            if not (read.has_tag(cb_tag) and read.has_tag(umi_tag)):
                continue
            cb = read.get_tag(cb_tag)
            ub = read.get_tag(umi_tag)
            if not (cb and ub):
                continue

            if require_gene:
                gene = get_gene_from_read(read, gene_tags)
                if gene is None:
                    if keep_missing_gene:
                        gene = "__NA__"
                    else:
                        continue
                item_key = (str(gene), str(ub))
            else:
                item_key = str(ub)

            cb = str(cb)
            cb2items[cb][item_key] += 1
            cb2N[cb] += 1

    return cb2items, cb2N


# ------------------------ sampling grids / combinatorics ----------------------

def make_sample_sizes(total_reads: int, min_reads: int, max_reads: int, num_points: int) -> np.ndarray:
    min_reads = max(2, min_reads)
    max_reads = min(max_reads, total_reads)
    if min_reads > max_reads:
        raise ValueError(f"min_reads ({min_reads}) > max_reads ({max_reads}); nothing to sample.")
    raw = np.logspace(math.log10(min_reads), math.log10(max_reads), num=num_points)
    sizes = np.unique(raw.astype(int))
    sizes = np.unique(np.concatenate([sizes, [min_reads, max_reads]])).astype(int)
    return sizes[(sizes >= min_reads) & (sizes <= max_reads)]


def make_sample_sizes_for_cell(N: int, global_min_reads: int, global_max_reads: Optional[int], num_points: int) -> np.ndarray:
    mr = max(2, global_min_reads)
    Mr = min(N, global_max_reads if global_max_reads is not None else N)
    if mr > Mr:
        mr, Mr = 2, N
    raw = np.logspace(math.log10(mr), math.log10(max(mr, Mr)), num=num_points)
    sizes = np.unique(raw.astype(int))
    sizes = np.unique(np.concatenate([sizes, [mr, Mr]])).astype(int)
    return sizes[(sizes >= 2) & (sizes <= N)]


# -------------------------- counting / expectations --------------------------

def count_unique(sampled: List[Tuple[str, str, Optional[str]]], dedup: str) -> int:
    if dedup == "cell":          # unique (CB, UMI)
        return len({(cb, ub) for cb, ub, _ in sampled})
    elif dedup == "library":     # unique UMI ignoring CB and Gene
        return len({ub for _, ub, _ in sampled})
    elif dedup == "cell-gene":   # unique (CB, Gene, UMI)
        return len({(cb, gene, ub) for cb, ub, gene in sampled})
    else:
        raise ValueError(f"Unknown --dedup: {dedup}")


def subsample_curve(records: List[Tuple[str, str, Optional[str]]],
                    sample_sizes: np.ndarray,
                    sample_times: int,
                    seed: int,
                    dedup: str) -> List[Tuple[int, float, float]]:
    random.seed(seed)
    n_total = len(records)
    results = []
    iterator = tqdm(sample_sizes, desc="Global subsampling", unit="pt") if HAS_TQDM else sample_sizes
    for n in iterator:
        uniq_counts = []
        for _ in range(sample_times):
            sampled = random.sample(records, n) if n < n_total else records
            uniq_counts.append(count_unique(sampled, dedup))
        mean_v = float(np.mean(uniq_counts))
        std_v = float(np.std(uniq_counts, ddof=1)) if len(uniq_counts) > 1 else 0.0
        results.append((int(n), mean_v, std_v))
    return results


def expected_unique_curve_for_cell(N: int,
                                   item_mult_counts: Counter,
                                   m_array: np.ndarray) -> np.ndarray:
    """
    Items are UMI (per-cell 'cell') or (Gene,UMI) (per-cell 'cell-gene').
    E = sum_c n_c * (1 - C(N-c, m) / C(N, m))
    """
    m = m_array.astype(int)
    logC_N_m = np.array([math.lgamma(N + 1) - math.lgamma(mm + 1) - math.lgamma(N - mm + 1) for mm in m])
    exp_unique = np.zeros_like(m, dtype=float)
    for c, n_c in item_mult_counts.items():
        term = np.ones_like(m, dtype=float)
        mask = m <= (N - c)
        if mask.any():
            logC_Nmc_m = np.array([
                math.lgamma(N - c + 1) - math.lgamma(mm + 1) - math.lgamma(N - c - mm + 1)
                for mm in m[mask]
            ])
            ratio = np.exp(logC_Nmc_m - logC_N_m[mask])
            term[mask] = 1.0 - ratio
        exp_unique += n_c * term
    return exp_unique


# --------------------------------- plotting ----------------------------------

def set_publication_style():
    mpl.rcParams.update({
        "figure.figsize": (7.5, 5.0),
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "axes.labelsize": 12,
        "axes.titlesize": 13,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
        "lines.linewidth": 2.0,
        "lines.markersize": 5,
        "font.family": "DejaVu Sans",
    })


def plot_global(results, out_png, title, xlabel, ylabel, legend_label, logx, show_band):
    set_publication_style()
    sizes = np.array([r[0] for r in results], dtype=float) / 1e6
    means = np.array([r[1] for r in results], dtype=float)
    stds  = np.array([r[2] for r in results], dtype=float)

    fig, ax = plt.subplots()
    ax.plot(sizes, means, marker="o", label=legend_label)
    if show_band:
        ax.fill_between(sizes, means - stds, means + stds, alpha=0.18, linewidth=0)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    if logx:
        ax.set_xscale("log")
        ax.get_xaxis().set_major_locator(mpl.ticker.LogLocator(base=10))
        ax.get_xaxis().set_minor_locator(mpl.ticker.LogLocator(base=10, subs=range(2, 10)))
        ax.get_xaxis().set_minor_formatter(mpl.ticker.NullFormatter())
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def plot_per_cell(curves, out_png, title, xlabel, ylabel, logx):
    set_publication_style()
    fig, ax = plt.subplots()
    for cb, (xs_reads, ys) in curves.items():
        ax.plot(xs_reads.astype(float) / 1e6, ys, linewidth=1.0, alpha=0.65)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    if logx:
        ax.set_xscale("log")
        ax.get_xaxis().set_major_locator(mpl.ticker.LogLocator(base=10))
        ax.get_xaxis().set_minor_locator(mpl.ticker.LogLocator(base=10, subs=range(2, 10)))
        ax.get_xaxis().set_minor_formatter(mpl.ticker.NullFormatter())
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


# ------------------------------- label helpers --------------------------------

def ylabel_default_for(mode: str, dedup: str) -> str:
    if mode == "global":
        if dedup == "cell":
            return "Unique (CB, UMI) count"
        elif dedup == "library":
            return "Unique UMI count"
        else:  # cell-gene
            return "Unique (CB, Gene, UMI) count"
    else:  # per-cell
        if dedup == "cell":
            return "Unique (CB, UMI) per cell (expected)"
        else:
            return "Unique (CB, Gene, UMI) per cell (expected)"


# ----------------------------------- CLI -------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Sequencing saturation curves for single-cell data."
    )
    p.add_argument("--mode", choices=["global", "per-cell"], default="global",
                   help="Which plot to generate (default: global).")

    # Common
    p.add_argument("--bam", required=True, help="Input BAM file.")
    p.add_argument("--cb-tag", default="CB", help="Cell barcode tag (default: CB).")
    p.add_argument("--umi-tag", default="UB", help="UMI tag (default: UB).")
    p.add_argument("--gene-tags", default="GX,GN",
                   help="Comma/space separated list of gene tags to try in order (default: GX,GN).")
    p.add_argument("--keep-missing-gene", action="store_true",
                   help="When a gene is required, keep reads without gene as '__NA__' instead of dropping.")
    p.add_argument("--min-mapq", type=int, default=0, help="Minimum MAPQ to keep reads (default: 0).")
    p.add_argument("--skip-secondary", action="store_true", help="Skip secondary/supplementary alignments.")
    p.add_argument("--logx", action="store_true", help="Use log-scale on the x-axis.")
    p.add_argument("--title", default=None, help="Custom plot title (optional).")
    p.add_argument("--xlabel", default=None, help="Override x-axis label.")
    p.add_argument("--ylabel", default=None, help="Override y-axis label.")
    p.add_argument("--out-prefix", required=True, help="Output prefix (PNG and optional CSV).")
    p.add_argument("--csv", action="store_true", help="Also write a CSV (long-format).")

    # Global
    p.add_argument("--dedup", choices=["cell", "library", "cell-gene"], default="cell",
                   help="Deduplication scope for global mode (default: cell).")
    p.add_argument("--num-points", type=int, default=25, help="Number of x-points (default: 25).")
    p.add_argument("--min-reads", type=int, default=5, help="Minimum reads to start subsampling (default: 5).")
    p.add_argument("--max-reads", type=int, default=None, help="Maximum reads to consider (default: use all).")
    p.add_argument("--sample-times", type=int, default=5, help="Monte Carlo repeats per point (default: 5).")
    p.add_argument("--no-band", action="store_true", help="Do not draw the ±SD band in global mode.")
    p.add_argument("--seed", type=int, default=42, help="Random seed (default: 42).")
    # Deprecated (compat)
    p.add_argument("--scope", choices=["per-cell", "across-cells"], help=argparse.SUPPRESS)

    # Per-cell
    p.add_argument("--per-cell-dedup", choices=["cell", "cell-gene"], default="cell",
                   help="Item definition for per-cell expected curves (default: cell).")
    p.add_argument("--per-cell-max", type=int, default=50,
                   help="Max number of cells to plot (top by reads; default: 50).")
    p.add_argument("--per-cell-min-reads", type=int, default=1000,
                   help="Only include cells with ≥ this many reads (default: 1000).")
    p.add_argument("--cb-list", type=str, default=None,
                   help="Optional file with CBs (one per line). If set, only these CBs are considered.")

    return p.parse_args()


def main():
    args = parse_args()

    # Backward compatibility for --scope
    if args.scope:
        print("[WARN] --scope is deprecated; use --dedup. Mapping automatically...")
        args.dedup = "cell" if args.scope == "per-cell" else "library"

    gene_tags = parse_gene_tags(args.gene_tags)
    default_xlabel = args.xlabel or "Number of reads (M)"

    # ---------------------------- PER-CELL MODE ----------------------------
    if args.mode == "per-cell":
        cb2items, cb2N = load_cb_item_counts(
            bam_path=args.bam,
            cb_tag=args.cb_tag,
            umi_tag=args.umi_tag,
            gene_tags=gene_tags,
            min_mapq=args.min_mapq,
            skip_secondary=args.skip_secondary,
            per_cell_dedup=args.per_cell_dedup,
            keep_missing_gene=args.keep_missing_gene,
        )

        # Filter & select cells
        if args.cb_list:
            with open(args.cb_list) as f:
                wanted = {ln.strip() for ln in f if ln.strip() and not ln.startswith("#")}
            selected = [(cb, cb2N[cb]) for cb in wanted if cb in cb2N and cb2N[cb] >= args.per_cell_min_reads]
        else:
            eligible = [(cb, n) for cb, n in cb2N.items() if n >= args.per_cell_min_reads]
            eligible.sort(key=lambda x: x[1], reverse=True)
            selected = eligible[:args.per_cell_max]

        if not selected:
            raise SystemExit("No cells selected. Check --cb-list / --per-cell-min-reads / data filters.")

        curves: Dict[str, Tuple[np.ndarray, np.ndarray]] = {}
        iterator = tqdm(selected, desc="Per-cell curves") if HAS_TQDM else selected
        for cb, N in iterator:
            mult_counter = Counter(cb2items[cb].values())
            if not mult_counter:
                continue
            m_array = make_sample_sizes_for_cell(
                N=N, global_min_reads=args.min_reads, global_max_reads=args.max_reads, num_points=args.num_points
            )
            ys = expected_unique_curve_for_cell(N, mult_counter, m_array)
            curves[cb] = (m_array, ys)

        if not curves:
            raise SystemExit("No per-cell curves to draw (no usable items). "
                             "Consider --keep-missing-gene or lowering thresholds.")

        out_png = f"{args.out_prefix}_per_cell_{args.per_cell_dedup}.png"
        ylabel_auto = ylabel_default_for("per-cell", args.per_cell_dedup)
        title = args.title or (f"Per-cell Saturation Curves (expected, {args.per_cell_dedup})\n"
                               f"{len(curves)} cells; min_reads≥{args.per_cell_min_reads}; "
                               f"{'CB list' if args.cb_list else f'top {args.per_cell_max} by reads'}")
        plot_per_cell(curves, out_png, title,
                      xlabel=default_xlabel,
                      ylabel=(args.ylabel or ylabel_auto),
                      logx=args.logx)
        print(f"[Per-cell] Figure: {out_png}")

        if args.csv:
            out_csv = f"{args.out_prefix}_per_cell_{args.per_cell_dedup}.csv"
            with open(out_csv, "w") as f:
                f.write("cell,reads,reads_M,expected_unique,dedup\n")
                for cb, (xs, ys) in curves.items():
                    for m, y in zip(xs, ys):
                        f.write(f"{cb},{int(m)},{m/1e6:.6f},{y:.6f},{args.per_cell_dedup}\n")
            print(f"[Per-cell] CSV:    {out_csv}")
        return

    # ----------------------------- GLOBAL MODE -----------------------------
    requires_gene = (args.dedup == "cell-gene")
    records = load_records(
        bam_path=args.bam,
        cb_tag=args.cb_tag,
        umi_tag=args.umi_tag,
        gene_tags=gene_tags,
        min_mapq=args.min_mapq,
        skip_secondary=args.skip_secondary,
        requires_gene=requires_gene,
        keep_missing_gene=args.keep_missing_gene,
    )
    total_reads = len(records)
    if total_reads == 0:
        msg = "No reads retained after filtering"
        if requires_gene and not args.keep_missing_gene:
            msg += " (try --keep-missing-gene or check --gene-tags)"
        raise SystemExit(msg + ".")

    max_reads = args.max_reads if args.max_reads is not None else total_reads
    sample_sizes = make_sample_sizes(total_reads, args.min_reads, max_reads, args.num_points)
    results = subsample_curve(records, sample_sizes, args.sample_times, seed=args.seed, dedup=args.dedup)

    if args.dedup == "cell":
        legend = "Unique (CB,UMI)"
        note = "(per-cell dedup)"
    elif args.dedup == "library":
        legend = "Unique UMI (library)"
        note = "(ignoring CB & Gene)"
    else:
        legend = "Unique (CB,Gene,UMI)"
        note = "(per-cell-per-gene dedup)"

    ylabel_auto = ylabel_default_for("global", args.dedup)
    out_png = f"{args.out_prefix}_global_{args.dedup}.png"
    title = args.title or f"Sequencing Saturation: Reads vs {legend}\nN={total_reads:,} reads {note}"
    plot_global(results, out_png, title,
                xlabel=default_xlabel,
                ylabel=(args.ylabel or ylabel_auto),
                legend_label=legend,
                logx=args.logx,
                show_band=not args.no_band)

    if args.csv:
        out_csv = f"{args.out_prefix}_global_{args.dedup}.csv"
        with open(out_csv, "w") as f:
            f.write("reads,reads_M,mean_unique,std_unique,dedup\n")
            for n, m, s in results:
                f.write(f"{n},{n/1e6:.6f},{m:.6f},{s:.6f},{args.dedup}\n")
        print(f"[Global]  CSV:    {out_csv}")

    max_pt = results[-1]
    print(f"[Global]  Parsed {total_reads:,} reads. "
          f"Peak: reads={max_pt[0]:,}, unique≈{int(round(max_pt[1])):,} ± {max_pt[2]:.2f} ({args.dedup}).")
    print(f"[Global]  Figure: {out_png}")


if __name__ == "__main__":
    main()
