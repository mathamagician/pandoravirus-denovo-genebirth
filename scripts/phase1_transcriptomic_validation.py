"""
Phase 1: Transcriptomic Validation of Poly-A/T Regulatory Signals
==================================================================
Pandoravirus 2.0

Compares Phase 9 poly-A/T regulatory predictions with published
transcriptomic data from Legendre et al. 2018.

Key questions (per peer review):
  1. Do poly-A/T runs coincide with experimental transcript boundaries?
  2. What is the relationship between hairpin terminators (Legendre 2018)
     and our poly-T terminator signals? Three hypotheses:
     a) Complementary: poly-A = promoters, hairpins = terminators
     b) Overlapping: poly-T + hairpin = rho-independent-like termination
     c) Our poly-T signal is a side-effect of AT-rich hairpin stems
  3. Does the 1.3-2.0 Mb "native zone" correlate with the 3' ORFan zone
     from the CRISPR study (Koonin/Abergel 2023)?
  4. Do proto-gene candidates (T4/T5) show any transcription?

Data sources:
  - Our database: genome sequence, gene positions, intergenic regions
  - Phase 9 code: poly-A/T run detection (reproduced here)
  - GBrowse reannotation GFF (if downloaded)
  - GenBank NC_022098.1: original annotation for comparison
"""

import sys
import os
import re
import time
import json
import numpy as np
import psycopg2
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from pathlib import Path
from collections import Counter, defaultdict
from scipy import stats as scipy_stats

# ── Paths ──────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parent.parent
ORIGINAL_PROJECT = Path("C:/Users/Eddie/Code Projects/Pandoravirus")
DATA_DIR = PROJECT_ROOT / "data"
FIGURES_DIR = PROJECT_ROOT / "figures"
RESULTS_DIR = PROJECT_ROOT / "results"

for d in [DATA_DIR, FIGURES_DIR, RESULTS_DIR]:
    d.mkdir(exist_ok=True)

DB_CONFIG = {
    "dbname": "pandoravirus",
    "user": "postgres",
    "password": "pandora2026",
    "host": "localhost",
    "port": 5432,
}


def get_conn():
    return psycopg2.connect(**DB_CONFIG)


# ═══════════════════════════════════════════════════════════════════════
# Data loading (matching Phase 9 patterns)
# ═══════════════════════════════════════════════════════════════════════

def get_sequence(cur, genome_id=4):
    """Retrieve full genome sequence from chunks."""
    cur.execute(
        "SELECT chunk_seq FROM sequence_chunk WHERE genome_id = %s ORDER BY chunk_start",
        (genome_id,))
    return ''.join(r[0] for r in cur.fetchall())


def get_genes(cur, genome_id=4):
    """Load genes with strand and product info."""
    cur.execute("""
        SELECT gene_id, start_pos, end_pos, strand, gene_length, product
        FROM gene WHERE genome_id = %s ORDER BY start_pos
    """, (genome_id,))
    cols = [d[0] for d in cur.description]
    return [dict(zip(cols, row)) for row in cur.fetchall()]


def is_hypothetical(product):
    if product is None:
        return True
    p = product.lower()
    return 'hypothetical' in p or p.startswith('orf') or 'unknown' in p


def find_poly_runs(seq, base, min_len=5):
    """Find all homopolymer runs of given base >= min_len."""
    runs = []
    pattern = re.compile(f'{base}{{{min_len},}}')
    for m in pattern.finditer(seq):
        runs.append({
            'start': m.start(),       # 0-based
            'end': m.end() - 1,       # 0-based inclusive
            'length': m.end() - m.start(),
            'pos_1based': m.start() + 1
        })
    return runs


COMPLEMENT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}

def revcomp(s):
    return ''.join(COMPLEMENT.get(b, 'N') for b in reversed(s))


# ═══════════════════════════════════════════════════════════════════════
# Analysis 1: Poly-A/T Inventory (reproduce Phase 9 baseline)
# ═══════════════════════════════════════════════════════════════════════

def poly_at_inventory(seq, genes):
    """Full characterization of poly-A and poly-T runs."""
    print("\n" + "=" * 70)
    print("ANALYSIS 1: Poly-A/T Run Inventory (Phase 9 Baseline)")
    print("=" * 70)

    # Build coding position set
    coding = set()
    for g in genes:
        for p in range(g['start_pos'], g['end_pos'] + 1):
            coding.add(p)

    results = {}
    for base in ['A', 'T']:
        runs = find_poly_runs(seq, base, min_len=5)
        results[base] = runs

        # Classify each run
        in_coding = 0
        in_intergenic = 0
        for r in runs:
            mid = r['pos_1based'] + r['length'] // 2
            if mid in coding:
                in_coding += 1
            else:
                in_intergenic += 1

        lengths = [r['length'] for r in runs]
        len_counts = Counter(lengths)
        ig_pct = in_intergenic / len(runs) * 100 if runs else 0

        print(f"\n  poly-{base} runs (>=5bp): {len(runs):>5,}")
        print(f"    Intergenic: {in_intergenic} ({ig_pct:.1f}%), Coding: {in_coding}")
        print(f"    Length: mean={np.mean(lengths):.1f}, median={np.median(lengths):.0f}, "
              f"max={max(lengths)}")
        print(f"    5bp={len_counts.get(5,0)}, 6bp={len_counts.get(6,0)}, "
              f"7bp={len_counts.get(7,0)}, 8+bp={sum(v for k,v in len_counts.items() if k>=8)}")

    return results


# ═══════════════════════════════════════════════════════════════════════
# Analysis 2: Hairpin Detection at Gene Boundaries
# ═══════════════════════════════════════════════════════════════════════

def find_hairpins(seq, position, window=50, min_stem=5, max_loop=10):
    """Search for palindromic/hairpin-forming sequences near a position.

    Legendre 2018: "70% of 3' transcript ends are overlapped by a
    nonconserved palindromic nucleotide sequence capable of forming a
    hairpin structure."

    Returns list of hairpins found within window of position.
    """
    start = max(0, position - window)
    end = min(len(seq), position + window)
    subseq = seq[start:end]

    hairpins = []
    for i in range(len(subseq) - 2 * min_stem):
        for stem_len in range(min_stem, min(15, (len(subseq) - i) // 2)):
            for loop_len in range(3, min(max_loop + 1, len(subseq) - i - 2 * stem_len + 1)):
                left = subseq[i:i + stem_len]
                right_start = i + stem_len + loop_len
                right_end = right_start + stem_len
                if right_end > len(subseq):
                    continue

                right = subseq[right_start:right_end]
                rc = revcomp(right)

                # Count complementary bases
                matches = sum(1 for a, b in zip(left, rc) if a == b)
                if matches >= stem_len * 0.8:  # Allow 20% mismatch
                    hairpins.append({
                        'genome_pos': start + i,
                        'stem_length': stem_len,
                        'loop_length': loop_len,
                        'total_length': 2 * stem_len + loop_len,
                        'distance_from_target': abs((start + i + stem_len) - position),
                        'matches': matches,
                        'mismatch_rate': 1 - matches / stem_len,
                    })

    # Keep only the best (longest stem) hairpin per position
    if hairpins:
        hairpins.sort(key=lambda h: (-h['stem_length'], h['distance_from_target']))
    return hairpins[:1] if hairpins else []


def hairpin_analysis(seq, genes):
    """Test relationship between hairpin structures and poly-T runs at gene 3' ends.

    Three hypotheses:
    a) Complementary: poly-A = promoters, hairpins = terminators (separate signals)
    b) Overlapping: poly-T + hairpin = rho-independent-like termination (co-localized)
    c) poly-T is a side-effect of AT-rich hairpin stems (confounded)
    """
    print("\n" + "=" * 70)
    print("ANALYSIS 2: Hairpin vs Poly-A/T Relationship at Gene 3' Ends")
    print("=" * 70)

    poly_T_runs = find_poly_runs(seq, 'T', min_len=5)

    # Build position lookup for poly-T runs
    poly_T_positions = set()
    for r in poly_T_runs:
        for p in range(r['start'], r['end'] + 1):
            poly_T_positions.add(p)

    # For each gene 3' end, check for hairpins and poly-T runs
    has_hairpin = 0
    has_polyT = 0
    has_both = 0
    has_neither = 0
    total = 0

    hairpin_only = []
    polyT_only = []
    both_signals = []

    for g in genes:
        # Gene 3' end (downstream of coding)
        if g['strand'] == '+':
            end_pos = g['end_pos']  # 1-based
        else:
            end_pos = g['start_pos']

        total += 1

        # Check for hairpin within 50bp downstream
        hp = find_hairpins(seq, end_pos - 1, window=50, min_stem=5)
        found_hp = len(hp) > 0

        # Check for poly-T within 100bp downstream
        if g['strand'] == '+':
            check_start = end_pos - 1
            check_end = min(len(seq), end_pos + 100)
        else:
            check_start = max(0, end_pos - 101)
            check_end = end_pos

        found_polyT = any(p in poly_T_positions for p in range(check_start, check_end))

        if found_hp and found_polyT:
            has_both += 1
            both_signals.append(g)
        elif found_hp:
            has_hairpin += 1
            hairpin_only.append(g)
        elif found_polyT:
            has_polyT += 1
            polyT_only.append(g)
        else:
            has_neither += 1

    print(f"\n  Gene 3' ends analyzed: {total}")
    print(f"    Hairpin + poly-T (both):  {has_both:>5} ({has_both/total*100:.1f}%)")
    print(f"    Hairpin only:             {has_hairpin:>5} ({has_hairpin/total*100:.1f}%)")
    print(f"    Poly-T only:              {has_polyT:>5} ({has_polyT/total*100:.1f}%)")
    print(f"    Neither:                  {has_neither:>5} ({has_neither/total*100:.1f}%)")

    # Interpretation
    total_with_hp = has_hairpin + has_both
    total_with_pT = has_polyT + has_both
    overlap_pct = has_both / max(1, min(total_with_hp, total_with_pT)) * 100

    print(f"\n  Total with hairpin:  {total_with_hp} ({total_with_hp/total*100:.1f}%)")
    print(f"  Total with poly-T:   {total_with_pT} ({total_with_pT/total*100:.1f}%)")
    print(f"  Overlap rate:        {overlap_pct:.1f}% of the smaller set")

    print(f"\n  --- Hypothesis Assessment ---")
    if overlap_pct > 60:
        print(f"  HYPOTHESIS B SUPPORTED: hairpin + poly-T co-localize")
        print(f"  --> rho-independent-like termination system")
    elif total_with_hp > total * 0.5 and total_with_pT < total * 0.3:
        print(f"  HYPOTHESIS A/C: hairpins dominate, poly-T is secondary")
    else:
        print(f"  MIXED: signals partially overlap, both contribute")

    return {
        'total': total,
        'hairpin_only': has_hairpin,
        'polyT_only': has_polyT,
        'both': has_both,
        'neither': has_neither,
        'total_hairpin': total_with_hp,
        'total_polyT': total_with_pT,
        'overlap_pct': overlap_pct,
    }


# ═══════════════════════════════════════════════════════════════════════
# Analysis 3: Spatial Distribution — 5'/3' Asymmetry (CRISPR correlation)
# ═══════════════════════════════════════════════════════════════════════

def spatial_orfan_analysis(genes, genome_length=2473870):
    """Test whether ORFan distribution correlates with CRISPR-identified
    3' ORFan zone from Koonin/Abergel 2023.

    Also check whether the 1.3-2.0 Mb "native zone" (Phase 2 GC dip)
    overlaps with the ORFan-enriched region.
    """
    print("\n" + "=" * 70)
    print("ANALYSIS 3: Spatial ORFan Distribution (CRISPR Correlation)")
    print("=" * 70)

    # Split genome into 10 equal bins
    n_bins = 10
    bin_size = genome_length / n_bins

    annotated_bins = np.zeros(n_bins)
    orfan_bins = np.zeros(n_bins)

    for g in genes:
        midpoint = (g['start_pos'] + g['end_pos']) / 2
        bin_idx = min(int(midpoint / bin_size), n_bins - 1)

        if is_hypothetical(g['product']):
            orfan_bins[bin_idx] += 1
        else:
            annotated_bins[bin_idx] += 1

    total_bins = annotated_bins + orfan_bins
    orfan_fraction = orfan_bins / np.maximum(total_bins, 1)

    print(f"\n  {'Bin':>6} {'Range (Mb)':>12} {'Annot':>6} {'ORFan':>6} {'Total':>6} {'ORFan%':>8}")
    print(f"  {'-'*50}")
    for i in range(n_bins):
        start_mb = i * bin_size / 1e6
        end_mb = (i + 1) * bin_size / 1e6
        print(f"  {i+1:>6} {start_mb:.1f}-{end_mb:.1f} Mb"
              f" {int(annotated_bins[i]):>6} {int(orfan_bins[i]):>6}"
              f" {int(total_bins[i]):>6} {orfan_fraction[i]*100:>7.1f}%")

    # Test 5' vs 3' enrichment
    # CRISPR paper: essential at 5', ORFan at 3'
    # Split at midpoint
    half = n_bins // 2
    orfan_5prime = orfan_bins[:half].sum()
    orfan_3prime = orfan_bins[half:].sum()
    total_5prime = total_bins[:half].sum()
    total_3prime = total_bins[half:].sum()

    frac_5 = orfan_5prime / total_5prime if total_5prime > 0 else 0
    frac_3 = orfan_3prime / total_3prime if total_3prime > 0 else 0

    print(f"\n  5' half ORFan fraction: {frac_5:.3f} ({int(orfan_5prime)}/{int(total_5prime)})")
    print(f"  3' half ORFan fraction: {frac_3:.3f} ({int(orfan_3prime)}/{int(total_3prime)})")

    # Chi-squared test
    contingency = np.array([
        [orfan_5prime, total_5prime - orfan_5prime],
        [orfan_3prime, total_3prime - orfan_3prime]
    ])
    chi2, p_val = scipy_stats.chi2_contingency(contingency)[:2]
    print(f"  Chi-squared: {chi2:.2f}, p = {p_val:.4f}")

    if p_val < 0.05 and frac_3 > frac_5:
        print(f"  --> CONSISTENT with CRISPR study: ORFans enriched at 3' end")
    elif p_val < 0.05 and frac_5 > frac_3:
        print(f"  --> OPPOSITE of CRISPR study: ORFans enriched at 5' end")
    else:
        print(f"  --> No significant 5'/3' asymmetry in ORFan distribution")

    # Check 1.3-2.0 Mb zone specifically
    zone_start_bin = int(1.3e6 / bin_size)
    zone_end_bin = min(int(2.0e6 / bin_size), n_bins - 1)
    zone_orfan = orfan_bins[zone_start_bin:zone_end_bin+1].sum()
    zone_total = total_bins[zone_start_bin:zone_end_bin+1].sum()
    zone_frac = zone_orfan / zone_total if zone_total > 0 else 0

    non_zone_orfan = orfan_bins.sum() - zone_orfan
    non_zone_total = total_bins.sum() - zone_total
    non_zone_frac = non_zone_orfan / non_zone_total if non_zone_total > 0 else 0

    print(f"\n  1.3-2.0 Mb zone: ORFan fraction = {zone_frac:.3f} ({int(zone_orfan)}/{int(zone_total)})")
    print(f"  Outside zone:    ORFan fraction = {non_zone_frac:.3f} ({int(non_zone_orfan)}/{int(non_zone_total)})")

    return {
        'annotated_bins': annotated_bins.tolist(),
        'orfan_bins': orfan_bins.tolist(),
        'orfan_fraction': orfan_fraction.tolist(),
        'frac_5prime': frac_5,
        'frac_3prime': frac_3,
        'chi2': chi2,
        'p_value': p_val,
        'zone_frac': zone_frac,
        'non_zone_frac': non_zone_frac,
    }


# ═══════════════════════════════════════════════════════════════════════
# Analysis 4: Gene Boundary Enrichment Test
# ═══════════════════════════════════════════════════════════════════════

def boundary_enrichment(seq, genes, poly_runs_A, poly_runs_T):
    """Enrichment of poly-A at gene starts vs poly-T at gene ends.

    This reproduces Phase 9's boundary analysis with additional
    distance distributions for comparison with transcriptomic data.
    """
    print("\n" + "=" * 70)
    print("ANALYSIS 4: Poly-A/T Boundary Enrichment (Extended)")
    print("=" * 70)

    # Strand-aware gene starts and ends
    gene_starts = []
    gene_ends = []
    for g in genes:
        if g['strand'] == '+':
            gene_starts.append(g['start_pos'])
            gene_ends.append(g['end_pos'])
        else:
            gene_starts.append(g['end_pos'])
            gene_ends.append(g['start_pos'])

    gene_starts = np.array(gene_starts)
    gene_ends = np.array(gene_ends)

    results = {}
    for label, runs, expected_near in [
        ('poly-A', poly_runs_A, 'starts'),
        ('poly-T', poly_runs_T, 'ends'),
    ]:
        positions = np.array([r['pos_1based'] for r in runs])

        # Distance to nearest gene start and end
        dist_to_starts = []
        dist_to_ends = []
        for pos in positions:
            d_start = np.min(np.abs(pos - gene_starts))
            d_end = np.min(np.abs(pos - gene_ends))
            dist_to_starts.append(d_start)
            dist_to_ends.append(d_end)

        dist_to_starts = np.array(dist_to_starts)
        dist_to_ends = np.array(dist_to_ends)

        # Count near boundaries (within 100bp)
        near_start_100 = np.sum(dist_to_starts <= 100)
        near_end_100 = np.sum(dist_to_ends <= 100)

        # Count near boundaries (within 50bp)
        near_start_50 = np.sum(dist_to_starts <= 50)
        near_end_50 = np.sum(dist_to_ends <= 50)

        print(f"\n  {label} runs (n={len(runs)}):")
        print(f"    Within 50bp of gene start:  {near_start_50:>5} ({near_start_50/len(runs)*100:.1f}%)")
        print(f"    Within 50bp of gene end:    {near_end_50:>5} ({near_end_50/len(runs)*100:.1f}%)")
        print(f"    Within 100bp of gene start: {near_start_100:>5} ({near_start_100/len(runs)*100:.1f}%)")
        print(f"    Within 100bp of gene end:   {near_end_100:>5} ({near_end_100/len(runs)*100:.1f}%)")
        print(f"    Median dist to start: {np.median(dist_to_starts):.0f} bp")
        print(f"    Median dist to end:   {np.median(dist_to_ends):.0f} bp")

        results[label] = {
            'n': len(runs),
            'near_start_50': int(near_start_50),
            'near_end_50': int(near_end_50),
            'near_start_100': int(near_start_100),
            'near_end_100': int(near_end_100),
            'median_dist_start': float(np.median(dist_to_starts)),
            'median_dist_end': float(np.median(dist_to_ends)),
            'dist_to_starts': dist_to_starts,
            'dist_to_ends': dist_to_ends,
        }

    # Test asymmetry: poly-A closer to starts? poly-T closer to ends?
    pA_closer_to_start = results['poly-A']['median_dist_start'] < results['poly-A']['median_dist_end']
    pT_closer_to_end = results['poly-T']['median_dist_end'] < results['poly-T']['median_dist_start']

    print(f"\n  --- Asymmetry Test ---")
    print(f"  poly-A closer to starts: {pA_closer_to_start} "
          f"(start={results['poly-A']['median_dist_start']:.0f} vs end={results['poly-A']['median_dist_end']:.0f})")
    print(f"  poly-T closer to ends:   {pT_closer_to_end} "
          f"(end={results['poly-T']['median_dist_end']:.0f} vs start={results['poly-T']['median_dist_start']:.0f})")

    if pA_closer_to_start and pT_closer_to_end:
        print(f"  --> CONFIRMED: Strand-asymmetric regulatory model")
        print(f"  --> poly-A = promoter signal, poly-T = terminator signal")
    else:
        print(f"  --> Asymmetry not confirmed at median level")

    return results


# ═══════════════════════════════════════════════════════════════════════
# Analysis 5: Rho-Independent Terminator Test
# ═══════════════════════════════════════════════════════════════════════

def rho_independent_test(seq, genes):
    """Test for rho-independent-like terminators (hairpin + poly-T tail).

    In bacteria, rho-independent termination uses a GC-rich hairpin
    followed by a poly-U (= poly-T in DNA) tail. If Pandoravirus uses
    a similar mechanism, we should find hairpins immediately followed
    by poly-T runs downstream of gene ends.
    """
    print("\n" + "=" * 70)
    print("ANALYSIS 5: Rho-Independent Terminator Structure Test")
    print("=" * 70)

    # For each gene, check 200bp downstream for hairpin+polyT pattern
    found_rho = 0
    found_hp_no_polyT = 0
    found_polyT_no_hp = 0
    found_neither = 0
    total = 0

    rho_details = []

    for g in genes:
        # Get downstream region (strand-aware)
        if g['strand'] == '+':
            ds_start = g['end_pos']  # 0-based = end_pos (since end_pos is 1-based inclusive)
            ds_end = min(len(seq), g['end_pos'] + 200)
        else:
            ds_start = max(0, g['start_pos'] - 201)
            ds_end = g['start_pos'] - 1

        if ds_end <= ds_start:
            continue

        downstream = seq[ds_start:ds_end]
        if g['strand'] == '-':
            downstream = revcomp(downstream)

        total += 1

        # Look for hairpin structure
        has_hairpin = False
        hp_end_pos = 0
        for i in range(min(100, len(downstream) - 15)):
            for stem in range(5, min(12, (len(downstream) - i) // 2)):
                for loop in range(3, 8):
                    if i + 2 * stem + loop > len(downstream):
                        continue
                    left = downstream[i:i + stem]
                    r_start = i + stem + loop
                    right = downstream[r_start:r_start + stem]
                    rc = revcomp(right)
                    matches = sum(1 for a, b in zip(left, rc) if a == b)
                    if matches >= stem * 0.8:
                        has_hairpin = True
                        hp_end_pos = r_start + stem
                        break
                if has_hairpin:
                    break
            if has_hairpin:
                break

        # Look for poly-T after hairpin (or anywhere in first 100bp)
        has_polyT = bool(re.search(r'T{5,}', downstream[:100]))

        # Specifically: poly-T within 20bp after hairpin end
        has_polyT_after_hp = False
        if has_hairpin and hp_end_pos < len(downstream):
            after_hp = downstream[hp_end_pos:min(hp_end_pos + 20, len(downstream))]
            has_polyT_after_hp = bool(re.search(r'T{4,}', after_hp))

        if has_hairpin and has_polyT_after_hp:
            found_rho += 1
            rho_details.append(g)
        elif has_hairpin:
            found_hp_no_polyT += 1
        elif has_polyT:
            found_polyT_no_hp += 1
        else:
            found_neither += 1

    print(f"\n  Genes analyzed: {total}")
    print(f"    Hairpin + poly-T tail (rho-like):  {found_rho:>5} ({found_rho/total*100:.1f}%)")
    print(f"    Hairpin only (no poly-T tail):     {found_hp_no_polyT:>5} ({found_hp_no_polyT/total*100:.1f}%)")
    print(f"    Poly-T only (no hairpin):          {found_polyT_no_hp:>5} ({found_polyT_no_hp/total*100:.1f}%)")
    print(f"    Neither:                           {found_neither:>5} ({found_neither/total*100:.1f}%)")

    total_with_hp = found_rho + found_hp_no_polyT
    print(f"\n  Total with hairpin: {total_with_hp} ({total_with_hp/total*100:.1f}%)")
    print(f"  (Legendre 2018 reported ~70% of 3' ends)")

    if found_rho > total * 0.3:
        print(f"\n  --> HYPOTHESIS B supported: substantial rho-independent-like termination")
        print(f"  --> Hairpin + poly-T co-occur at {found_rho/total*100:.1f}% of gene ends")
    elif total_with_hp > total * 0.5 and found_polyT_no_hp > total * 0.15:
        print(f"\n  --> MIXED: Both mechanisms present, partially independent")
    else:
        print(f"\n  --> Pattern unclear, further analysis needed")

    return {
        'total': total,
        'rho_like': found_rho,
        'hairpin_only': found_hp_no_polyT,
        'polyT_only': found_polyT_no_hp,
        'neither': found_neither,
    }


# ═══════════════════════════════════════════════════════════════════════
# Visualization
# ═══════════════════════════════════════════════════════════════════════

def plot_phase1_results(hairpin_data, spatial_data, boundary_data, rho_data):
    """4-panel Phase 1 figure."""
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Phase 1: Transcriptomic Validation -- Pandoravirus 2.0',
                 fontsize=14, fontweight='bold')

    # Panel A: Hairpin vs Poly-T at gene 3' ends
    ax = axes[0, 0]
    labels = ['Hairpin +\npoly-T', 'Hairpin\nonly', 'Poly-T\nonly', 'Neither']
    sizes = [hairpin_data['both'], hairpin_data['hairpin_only'],
             hairpin_data['polyT_only'], hairpin_data['neither']]
    colors = ['#2ecc71', '#3498db', '#e74c3c', '#95a5a6']
    bars = ax.bar(labels, sizes, color=colors, edgecolor='black', alpha=0.8)
    ax.set_ylabel('Number of gene 3\' ends')
    ax.set_title('A. Hairpin vs Poly-T at Gene 3\' Ends')
    for bar, size in zip(bars, sizes):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 5,
                f'{size}\n({size/hairpin_data["total"]*100:.0f}%)',
                ha='center', va='bottom', fontsize=8)

    # Panel B: Spatial ORFan distribution
    ax = axes[0, 1]
    n_bins = len(spatial_data['orfan_fraction'])
    x = np.arange(n_bins)
    width = 0.35
    ax.bar(x - width/2, spatial_data['annotated_bins'], width,
           label='Annotated', color='steelblue', alpha=0.8)
    ax.bar(x + width/2, spatial_data['orfan_bins'], width,
           label='ORFan', color='coral', alpha=0.8)
    ax.axvspan(5.2, 8.0, alpha=0.1, color='green', label='1.3-2.0 Mb zone')
    ax.set_xlabel('Genome position (deciles)')
    ax.set_ylabel('Gene count')
    ax.set_title('B. Spatial Gene Distribution (5\'->3\')')
    ax.legend(fontsize=8)
    bin_labels = [f'{i*0.247:.1f}' for i in range(n_bins)]
    ax.set_xticks(x)
    ax.set_xticklabels(bin_labels, fontsize=7)
    ax.text(0.02, 0.95, f'p={spatial_data["p_value"]:.4f}',
            transform=ax.transAxes, fontsize=8, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

    # Panel C: Distance distributions
    ax = axes[1, 0]
    if 'poly-A' in boundary_data:
        dist_start = boundary_data['poly-A']['dist_to_starts']
        dist_end = boundary_data['poly-T']['dist_to_ends']
        ax.hist(dist_start[dist_start <= 500], bins=50, alpha=0.6, color='royalblue',
                label=f'poly-A to starts (med={np.median(dist_start):.0f})')
        ax.hist(dist_end[dist_end <= 500], bins=50, alpha=0.6, color='tomato',
                label=f'poly-T to ends (med={np.median(dist_end):.0f})')
    ax.set_xlabel('Distance to nearest boundary (bp)')
    ax.set_ylabel('Count')
    ax.set_title('C. Poly-A/T Distance to Gene Boundaries')
    ax.legend(fontsize=8)

    # Panel D: Rho-independent terminator structure
    ax = axes[1, 1]
    labels = ['Hairpin +\npoly-T tail\n(rho-like)',
              'Hairpin\nonly', 'Poly-T\nonly', 'Neither']
    sizes = [rho_data['rho_like'], rho_data['hairpin_only'],
             rho_data['polyT_only'], rho_data['neither']]
    colors = ['#2ecc71', '#3498db', '#e74c3c', '#95a5a6']
    ax.pie(sizes, labels=labels, autopct='%1.1f%%', colors=colors,
           textprops={'fontsize': 8})
    ax.set_title('D. Termination Mechanism Classification')

    plt.tight_layout()
    outpath = FIGURES_DIR / 'phase1_transcriptomic_validation.png'
    plt.savefig(str(outpath), dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\n  Figure saved: {outpath}")


# ═══════════════════════════════════════════════════════════════════════
# MAIN
# ═══════════════════════════════════════════════════════════════════════

def main():
    t0 = time.time()

    print()
    print("=" * 70)
    print("PANDORAVIRUS 2.0 -- PHASE 1: TRANSCRIPTOMIC VALIDATION")
    print("=" * 70)

    # Connect to database
    conn = get_conn()
    cur = conn.cursor()

    # Load genome sequence
    print("\nLoading Pandoravirus salinus genome...")
    seq = get_sequence(cur, genome_id=4)
    print(f"  Genome length: {len(seq):,} bp")

    # Load genes
    genes = get_genes(cur, genome_id=4)
    print(f"  Genes: {len(genes)}")

    annotated = [g for g in genes if not is_hypothetical(g['product'])]
    orfan = [g for g in genes if is_hypothetical(g['product'])]
    print(f"  Annotated: {len(annotated)}, ORFan/Hypothetical: {len(orfan)}")

    # Analysis 1: Poly-A/T inventory (baseline)
    poly_runs = poly_at_inventory(seq, genes)

    # Analysis 2: Hairpin vs poly-A/T at 3' ends
    hairpin_data = hairpin_analysis(seq, genes)

    # Analysis 3: Spatial ORFan distribution (CRISPR correlation)
    spatial_data = spatial_orfan_analysis(genes)

    # Analysis 4: Boundary enrichment with distance distributions
    boundary_data = boundary_enrichment(seq, genes, poly_runs['A'], poly_runs['T'])

    # Analysis 5: Rho-independent terminator test
    rho_data = rho_independent_test(seq, genes)

    # Visualization
    print("\n" + "=" * 70)
    print("GENERATING FIGURES")
    print("=" * 70)
    # Remove numpy arrays before saving (not JSON serializable)
    boundary_save = {}
    for key in boundary_data:
        boundary_save[key] = {k: v for k, v in boundary_data[key].items()
                              if k not in ('dist_to_starts', 'dist_to_ends')}

    plot_phase1_results(hairpin_data, spatial_data, boundary_data, rho_data)

    # Save results
    results = {
        'hairpin_analysis': hairpin_data,
        'spatial_distribution': spatial_data,
        'boundary_enrichment': boundary_save,
        'rho_independent': rho_data,
    }

    results_path = RESULTS_DIR / 'phase1_results.json'
    with open(str(results_path), 'w') as f:
        json.dump(results, f, indent=2)
    print(f"  Results saved: {results_path}")

    # Summary
    elapsed = time.time() - t0
    print(f"\n" + "=" * 70)
    print(f"PHASE 1 SUMMARY")
    print(f"=" * 70)
    print(f"  Runtime: {elapsed:.1f}s")
    print(f"\n  Key findings:")
    print(f"  1. Hairpin + poly-T co-occurrence at 3' ends: "
          f"{hairpin_data['both']}/{hairpin_data['total']} ({hairpin_data['both']/hairpin_data['total']*100:.1f}%)")
    print(f"  2. 5'/3' ORFan asymmetry: p={spatial_data['p_value']:.4f}")
    print(f"     5' ORFan fraction: {spatial_data['frac_5prime']:.3f}")
    print(f"     3' ORFan fraction: {spatial_data['frac_3prime']:.3f}")
    print(f"  3. 1.3-2.0 Mb zone ORFan fraction: {spatial_data['zone_frac']:.3f} "
          f"vs outside: {spatial_data['non_zone_frac']:.3f}")
    print(f"  4. Poly-A/T boundary asymmetry confirmed: "
          f"{'YES' if boundary_data['poly-A']['median_dist_start'] < boundary_data['poly-A']['median_dist_end'] else 'NO'}")

    print(f"\n  Next steps:")
    print(f"  - Download reannotated GFF from GBrowse for transcript boundary comparison")
    print(f"  - If hairpin+polyT co-localize strongly: rho-independent model")
    print(f"  - Proceed to Phase 2 (BLASTp) in parallel")

    conn.close()


if __name__ == "__main__":
    main()
