"""
Phase 1b: Corrected Analysis with Three Fixes
==============================================
Pandoravirus 2.0

Fixes applied:
  1. Intergenic-only filter on poly-A/T runs (preprocessing)
  2. Strand-aware boundary test (reproduces Phase 9 meta-gene asymmetry)
  3. Hairpin null model (Markov-1 shuffled sequence calibration)
  4. Two-population test (active vs latent regulatory signals)
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
FIGURES_DIR = PROJECT_ROOT / "figures"
RESULTS_DIR = PROJECT_ROOT / "results"

for d in [FIGURES_DIR, RESULTS_DIR]:
    d.mkdir(exist_ok=True)

DB_CONFIG = {
    "dbname": "pandoravirus",
    "user": "postgres",
    "password": "pandora2026",
    "host": "localhost",
    "port": 5432,
}

COMPLEMENT = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}

def revcomp(s):
    return ''.join(COMPLEMENT.get(b, 'N') for b in reversed(s))


def get_conn():
    return psycopg2.connect(**DB_CONFIG)


def get_sequence(cur, genome_id=4):
    cur.execute(
        "SELECT chunk_seq FROM sequence_chunk WHERE genome_id = %s ORDER BY chunk_start",
        (genome_id,))
    return ''.join(r[0] for r in cur.fetchall())


def get_genes(cur, genome_id=4):
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
    runs = []
    pattern = re.compile(f'{base}{{{min_len},}}')
    for m in pattern.finditer(seq):
        runs.append({
            'start': m.start(),
            'end': m.end() - 1,
            'length': m.end() - m.start(),
            'pos_1based': m.start() + 1,
            'midpoint': (m.start() + m.end()) // 2,
        })
    return runs


# ═══════════════════════════════════════════════════════════════════════
# FIX 1: Intergenic-Only Filter
# ═══════════════════════════════════════════════════════════════════════

def filter_intergenic_runs(runs, genes):
    """Keep only poly-A/T runs whose midpoint falls in intergenic space."""
    coding = set()
    for g in genes:
        for p in range(g['start_pos'], g['end_pos'] + 1):
            coding.add(p)

    intergenic_runs = []
    coding_runs = []
    for r in runs:
        mid = r['pos_1based'] + r['length'] // 2
        if mid in coding:
            coding_runs.append(r)
        else:
            intergenic_runs.append(r)

    return intergenic_runs, coding_runs


# ═══════════════════════════════════════════════════════════════════════
# FIX 2: Strand-Aware Boundary Test
# ═══════════════════════════════════════════════════════════════════════

def strand_aware_boundary_analysis(seq, genes, ig_poly_A, ig_poly_T, window=200):
    """Reproduce Phase 9's strand-aware asymmetry test.

    For + strand genes:
      - Upstream = lower coordinates (before start_pos)
      - Downstream = higher coordinates (after end_pos)
      - Poly-A upstream of start = promoter signal
      - Poly-T downstream of end = terminator signal

    For - strand genes:
      - Upstream = higher coordinates (after end_pos)
      - Downstream = lower coordinates (before start_pos)
      - Poly-A upstream of end = promoter signal (on reverse strand)
      - Poly-T downstream of start = terminator signal (on reverse strand)

    Actually, on the - strand, poly-T on the genome = poly-A on template.
    So for - strand: poly-T upstream = poly-A on template = promoter.
    This means strand-aware analysis needs to account for this.
    """
    print("\n" + "=" * 70)
    print("FIX 2: Strand-Aware Boundary Test")
    print("=" * 70)

    # Gene start = transcription start (strand-aware)
    # Gene end = transcription end (strand-aware)
    plus_genes = [g for g in genes if g['strand'] == '+']
    minus_genes = [g for g in genes if g['strand'] == '-']

    poly_A_pos = np.array([r['pos_1based'] for r in ig_poly_A])
    poly_T_pos = np.array([r['pos_1based'] for r in ig_poly_T])

    results = {}

    # ── Plus strand analysis ──
    # Promoter test: poly-A runs upstream of + strand gene starts
    plus_starts = np.array([g['start_pos'] for g in plus_genes])
    plus_ends = np.array([g['end_pos'] for g in plus_genes])

    # Count poly-A within window bp upstream of + strand starts
    pA_upstream_of_plus_start = 0
    pA_downstream_of_plus_start = 0
    pT_upstream_of_plus_end = 0
    pT_downstream_of_plus_end = 0

    for pos in poly_A_pos:
        for gs in plus_starts:
            if gs - window <= pos < gs:   # upstream
                pA_upstream_of_plus_start += 1
                break
            elif gs <= pos < gs + window:  # downstream (inside gene)
                pA_downstream_of_plus_start += 1
                break

    for pos in poly_T_pos:
        for ge in plus_ends:
            if ge - window < pos <= ge:    # upstream (inside gene)
                pT_upstream_of_plus_end += 1
                break
            elif ge < pos <= ge + window:  # downstream
                pT_downstream_of_plus_end += 1
                break

    # ── Minus strand analysis ──
    # For - strand: transcription starts at end_pos, ends at start_pos
    # On the template (- strand), poly-T on genome = poly-A on template
    # So poly-T upstream of - strand start (= high coord) = promoter equivalent
    minus_starts = np.array([g['end_pos'] for g in minus_genes])  # transcription start
    minus_ends = np.array([g['start_pos'] for g in minus_genes])  # transcription end

    pT_upstream_of_minus_start = 0   # poly-T above end_pos = poly-A on template = promoter
    pT_downstream_of_minus_start = 0
    pA_upstream_of_minus_end = 0     # poly-A below start_pos = poly-T on template = terminator
    pA_downstream_of_minus_end = 0

    for pos in poly_T_pos:
        for gs in minus_starts:
            if gs < pos <= gs + window:   # upstream on - strand (higher coords)
                pT_upstream_of_minus_start += 1
                break
            elif gs - window <= pos <= gs:  # downstream (inside gene)
                pT_downstream_of_minus_start += 1
                break

    for pos in poly_A_pos:
        for ge in minus_ends:
            if ge - window <= pos < ge:    # upstream on - strand (lower coords = inside gene)
                pA_upstream_of_minus_end += 1
                break
            elif ge - window <= pos < ge:
                pA_downstream_of_minus_end += 1
                break

    # Combined: promoter = poly-A upstream of + starts + poly-T upstream of - starts
    # Combined: terminator = poly-T downstream of + ends + poly-A downstream of - ends
    promoter_signal = pA_upstream_of_plus_start + pT_upstream_of_minus_start
    terminator_signal = pT_downstream_of_plus_end + pA_downstream_of_minus_end

    # Control: wrong-signal combinations
    wrong_promoter = pT_upstream_of_plus_end + pA_upstream_of_minus_end  # poly-T before start? Wrong.
    wrong_terminator = pA_downstream_of_plus_start + pT_downstream_of_minus_start  # poly-A after end? Wrong.

    print(f"\n  Plus strand genes: {len(plus_genes)}")
    print(f"  Minus strand genes: {len(minus_genes)}")
    print(f"  Intergenic poly-A runs: {len(ig_poly_A)}")
    print(f"  Intergenic poly-T runs: {len(ig_poly_T)}")

    print(f"\n  --- Strand-Aware Boundary Counts (window={window}bp) ---")
    print(f"  Plus strand:")
    print(f"    poly-A upstream of start (promoter):     {pA_upstream_of_plus_start}")
    print(f"    poly-A downstream of start (into gene):  {pA_downstream_of_plus_start}")
    print(f"    poly-T downstream of end (terminator):   {pT_downstream_of_plus_end}")
    print(f"    poly-T upstream of end (into gene):      {pT_upstream_of_plus_end}")

    print(f"  Minus strand:")
    print(f"    poly-T upstream of start (promoter):     {pT_upstream_of_minus_start}")
    print(f"    poly-T downstream of start (into gene):  {pT_downstream_of_minus_start}")
    print(f"    poly-A downstream of end (terminator):   {pA_downstream_of_minus_end}")

    print(f"\n  --- Combined (strand-corrected) ---")
    print(f"    PROMOTER signal (correct):   {promoter_signal}")
    print(f"    TERMINATOR signal (correct): {terminator_signal}")

    # Chi-squared test: correct signal vs wrong position
    total_correct = promoter_signal + terminator_signal
    total_all_boundary = total_correct + wrong_promoter + wrong_terminator

    print(f"\n    Total boundary-associated (correct): {total_correct}")
    print(f"    Total boundary-associated (all):     {total_all_boundary}")

    if total_all_boundary > 0:
        correct_frac = total_correct / total_all_boundary
        print(f"    Correct fraction: {correct_frac:.3f}")

        # Binomial test: is correct fraction > 0.5?
        from scipy.stats import binomtest
        binom_result = binomtest(total_correct, total_all_boundary, 0.5, alternative='greater')
        print(f"    Binomial test (correct > 50%): p = {binom_result.pvalue:.6f}")

    results = {
        'plus_genes': len(plus_genes),
        'minus_genes': len(minus_genes),
        'pA_upstream_plus_start': pA_upstream_of_plus_start,
        'pT_downstream_plus_end': pT_downstream_of_plus_end,
        'pT_upstream_minus_start': pT_upstream_of_minus_start,
        'promoter_signal': promoter_signal,
        'terminator_signal': terminator_signal,
        'total_correct': total_correct,
        'total_all_boundary': total_all_boundary,
    }

    return results


def metagene_at_profile(seq, genes, window=500, step=10):
    """Compute average AT-content profile around gene starts and ends.

    Exactly reproduces Phase 9's meta-gene AT profile.
    Positive offset = into gene, negative = into intergenic.
    """
    seq_len = len(seq)
    positions = np.arange(-window, window + 1, step)

    start_profiles = []
    end_profiles = []

    for g in genes:
        # Gene start (TSS proxy)
        if g['strand'] == '+':
            anchor = g['start_pos'] - 1  # 0-based
        else:
            anchor = g['end_pos'] - 1

        profile = []
        for offset in positions:
            if g['strand'] == '+':
                pos = anchor + offset
            else:
                pos = anchor - offset
            if 0 <= pos < seq_len:
                profile.append(1.0 if seq[pos] in 'AT' else 0.0)
            else:
                profile.append(np.nan)
        start_profiles.append(profile)

        # Gene end
        if g['strand'] == '+':
            anchor_end = g['end_pos'] - 1
        else:
            anchor_end = g['start_pos'] - 1

        profile_end = []
        for offset in positions:
            if g['strand'] == '+':
                pos = anchor_end + offset
            else:
                pos = anchor_end - offset
            if 0 <= pos < seq_len:
                profile_end.append(1.0 if seq[pos] in 'AT' else 0.0)
            else:
                profile_end.append(np.nan)
        end_profiles.append(profile_end)

    start_arr = np.nanmean(start_profiles, axis=0)
    end_arr = np.nanmean(end_profiles, axis=0)

    return positions, start_arr, end_arr


# ═══════════════════════════════════════════════════════════════════════
# FIX 3: Hairpin Null Model
# ═══════════════════════════════════════════════════════════════════════

def generate_markov1_sequence(seq, rng):
    """Generate a Markov-1 (dinucleotide-preserving) shuffled sequence."""
    # Count dinucleotide transitions
    transitions = defaultdict(lambda: defaultdict(int))
    for i in range(len(seq) - 1):
        transitions[seq[i]][seq[i+1]] += 1

    # Convert to probabilities
    trans_prob = {}
    for base, nexts in transitions.items():
        total = sum(nexts.values())
        trans_prob[base] = {n: c/total for n, c in nexts.items()}

    # Generate sequence
    result = [seq[0]]  # Start with same first base
    for _ in range(len(seq) - 1):
        current = result[-1]
        if current in trans_prob:
            bases = list(trans_prob[current].keys())
            probs = list(trans_prob[current].values())
            next_base = rng.choice(bases, p=probs)
            result.append(next_base)
        else:
            result.append(rng.choice(['A', 'C', 'G', 'T']))

    return ''.join(result)


def find_hairpins_in_region(subseq, min_stem=5, max_loop=10, max_stem=15):
    """Find the best hairpin in a subsequence.

    Returns True if a hairpin with >= min_stem complementary bases found.
    Uses stricter matching than Phase 1a (no mismatch tolerance).
    """
    for i in range(len(subseq) - 2 * min_stem):
        for stem_len in range(min(max_stem, (len(subseq) - i) // 2), min_stem - 1, -1):
            for loop_len in range(3, min(max_loop + 1, len(subseq) - i - 2 * stem_len + 1)):
                left = subseq[i:i + stem_len]
                r_start = i + stem_len + loop_len
                right = subseq[r_start:r_start + stem_len]
                rc = revcomp(right)

                # STRICT: require exact complementarity (no mismatches)
                if left == rc:
                    return True
    return False


def hairpin_null_model(seq, genes, n_shuffles=20, min_stem=6):
    """Calibrate hairpin detection against Markov-1 null.

    Tests whether hairpin frequency at gene 3' ends exceeds
    what's expected from sequence composition alone.

    Uses stricter parameters than Phase 1a:
      - min_stem=6 (was 5)
      - Exact complementarity (was 80% match)
    """
    print("\n" + "=" * 70)
    print("FIX 3: Hairpin Null Model (Markov-1 Calibration)")
    print("=" * 70)
    print(f"  Parameters: min_stem={min_stem}, exact complementarity, window=60bp")
    print(f"  Null shuffles: {n_shuffles}")

    rng = np.random.default_rng(42)

    # Extract 60bp downstream of each gene 3' end
    regions = []
    for g in genes:
        if g['strand'] == '+':
            ds_start = g['end_pos']
            ds_end = min(len(seq), g['end_pos'] + 60)
        else:
            ds_start = max(0, g['start_pos'] - 61)
            ds_end = g['start_pos'] - 1

        if ds_end - ds_start < 15:
            continue

        subseq = seq[ds_start:ds_end]
        if g['strand'] == '-':
            subseq = revcomp(subseq)
        regions.append(subseq)

    # Observed hairpin rate
    observed_count = sum(1 for r in regions if find_hairpins_in_region(r, min_stem=min_stem))
    observed_rate = observed_count / len(regions)
    print(f"\n  Observed hairpin rate: {observed_count}/{len(regions)} = {observed_rate:.3f} ({observed_rate*100:.1f}%)")

    # Null: shuffle each region (Markov-1) and count hairpins
    null_rates = []
    for i in range(n_shuffles):
        null_count = 0
        for r in regions:
            shuffled = generate_markov1_sequence(r, rng)
            if find_hairpins_in_region(shuffled, min_stem=min_stem):
                null_count += 1
        null_rate = null_count / len(regions)
        null_rates.append(null_rate)
        if (i + 1) % 5 == 0:
            print(f"    Shuffle {i+1}/{n_shuffles}: null rate = {null_rate:.3f}")

    null_rates = np.array(null_rates)
    null_mean = null_rates.mean()
    null_std = null_rates.std()
    z_score = (observed_rate - null_mean) / null_std if null_std > 0 else 0
    p_value = np.mean(null_rates >= observed_rate)

    print(f"\n  Null mean: {null_mean:.3f} ({null_mean*100:.1f}%)")
    print(f"  Null std:  {null_std:.3f}")
    print(f"  Z-score:   {z_score:.2f}")
    print(f"  p-value:   {p_value:.4f}")

    excess = observed_rate - null_mean
    print(f"\n  Excess over null: {excess:.3f} ({excess*100:.1f} percentage points)")
    print(f"  Legendre 2018 reported ~70% -- our strict rate: {observed_rate*100:.1f}%")

    if z_score > 2:
        print(f"  --> SIGNIFICANT: Hairpins enriched beyond sequence composition (z={z_score:.1f})")
        print(f"  --> But excess is {excess*100:.1f}pp -- {('substantial' if excess > 0.1 else 'modest')} biological signal")
    else:
        print(f"  --> NOT SIGNIFICANT: Hairpin rate explained by sequence composition")
        print(f"  --> Our Phase 1a 99.5% detection was indeed noise")

    return {
        'observed_rate': observed_rate,
        'observed_count': observed_count,
        'n_regions': len(regions),
        'null_mean': null_mean,
        'null_std': null_std,
        'z_score': z_score,
        'p_value': p_value,
        'excess': excess,
        'min_stem': min_stem,
    }


# ═══════════════════════════════════════════════════════════════════════
# Analysis: Two-Population Test (Active vs Latent Regulatory Signals)
# ═══════════════════════════════════════════════════════════════════════

def two_population_test(ig_poly_A, ig_poly_T, genes, boundary_dist=100):
    """Test whether boundary-proximal poly-A/T runs differ from dispersed ones.

    Hypothesis: boundary-proximal runs are "active" regulators (should show
    strand asymmetry). Dispersed runs are "latent" regulators (no strand
    preference, but provide regulatory-ready substrate for gene birth).
    """
    print("\n" + "=" * 70)
    print("ANALYSIS: Two-Population Test (Active vs Latent Signals)")
    print("=" * 70)

    # Gene boundary positions (strand-aware starts and ends)
    starts = []
    ends = []
    for g in genes:
        if g['strand'] == '+':
            starts.append(g['start_pos'])
            ends.append(g['end_pos'])
        else:
            starts.append(g['end_pos'])
            ends.append(g['start_pos'])

    starts = np.array(starts)
    ends = np.array(ends)

    def classify_run(run_pos):
        """Classify a run as boundary-proximal or dispersed."""
        min_dist_start = np.min(np.abs(run_pos - starts))
        min_dist_end = np.min(np.abs(run_pos - ends))
        min_dist = min(min_dist_start, min_dist_end)

        if min_dist <= boundary_dist:
            return 'boundary', min_dist, min_dist_start, min_dist_end
        else:
            return 'dispersed', min_dist, min_dist_start, min_dist_end

    # Classify all intergenic poly-A and poly-T runs
    boundary_A = []
    dispersed_A = []
    boundary_T = []
    dispersed_T = []

    for r in ig_poly_A:
        cat, dist, ds, de = classify_run(r['pos_1based'])
        if cat == 'boundary':
            boundary_A.append(r)
        else:
            dispersed_A.append(r)

    for r in ig_poly_T:
        cat, dist, ds, de = classify_run(r['pos_1based'])
        if cat == 'boundary':
            boundary_T.append(r)
        else:
            dispersed_T.append(r)

    print(f"\n  Boundary-proximal (<={boundary_dist}bp from gene boundary):")
    print(f"    poly-A: {len(boundary_A):>5}  poly-T: {len(boundary_T):>5}")
    print(f"  Dispersed (>{boundary_dist}bp from any boundary):")
    print(f"    poly-A: {len(dispersed_A):>5}  poly-T: {len(dispersed_T):>5}")

    # Test 1: Length distribution difference
    boundary_lengths = [r['length'] for r in boundary_A + boundary_T]
    dispersed_lengths = [r['length'] for r in dispersed_A + dispersed_T]

    if boundary_lengths and dispersed_lengths:
        u, p = scipy_stats.mannwhitneyu(boundary_lengths, dispersed_lengths, alternative='two-sided')
        print(f"\n  Length comparison:")
        print(f"    Boundary mean: {np.mean(boundary_lengths):.2f} bp")
        print(f"    Dispersed mean: {np.mean(dispersed_lengths):.2f} bp")
        print(f"    Mann-Whitney p = {p:.4f}")

    # Test 2: A/T ratio in boundary vs dispersed
    boundary_A_frac = len(boundary_A) / (len(boundary_A) + len(boundary_T)) if (len(boundary_A) + len(boundary_T)) > 0 else 0
    dispersed_A_frac = len(dispersed_A) / (len(dispersed_A) + len(dispersed_T)) if (len(dispersed_A) + len(dispersed_T)) > 0 else 0

    print(f"\n  A/T ratio:")
    print(f"    Boundary poly-A fraction: {boundary_A_frac:.3f}")
    print(f"    Dispersed poly-A fraction: {dispersed_A_frac:.3f}")

    # Chi-squared test on A vs T proportion
    contingency = np.array([
        [len(boundary_A), len(boundary_T)],
        [len(dispersed_A), len(dispersed_T)]
    ])
    chi2, p_chi = scipy_stats.chi2_contingency(contingency)[:2]
    print(f"    Chi-squared: {chi2:.2f}, p = {p_chi:.4f}")

    if p_chi < 0.05:
        print(f"    --> Boundary and dispersed populations have DIFFERENT A/T compositions")
    else:
        print(f"    --> A/T composition is SIMILAR between populations")

    return {
        'n_boundary_A': len(boundary_A),
        'n_boundary_T': len(boundary_T),
        'n_dispersed_A': len(dispersed_A),
        'n_dispersed_T': len(dispersed_T),
        'boundary_A_frac': boundary_A_frac,
        'dispersed_A_frac': dispersed_A_frac,
        'length_p': p if boundary_lengths and dispersed_lengths else None,
        'AT_ratio_p': p_chi,
    }


# ═══════════════════════════════════════════════════════════════════════
# Visualization
# ═══════════════════════════════════════════════════════════════════════

def plot_phase1b_results(positions, start_profile, end_profile,
                         strand_data, hairpin_data, twopop_data,
                         ig_poly_A=None, ig_poly_T=None, genes=None):
    """6-panel corrected figure."""
    fig, axes = plt.subplots(2, 3, figsize=(18, 10))
    fig.suptitle('Phase 1b: Corrected Transcriptomic Validation -- Pandoravirus 2.0',
                 fontsize=14, fontweight='bold')

    # Panel A: Meta-gene AT profile (reproduced from Phase 9)
    ax = axes[0, 0]
    genome_at = (1 - 0.6172) * 100  # 38.28% AT (from GC=61.72%)
    ax.plot(positions, start_profile * 100, 'b-', linewidth=1.5, label='Gene START', alpha=0.9)
    ax.plot(positions, end_profile * 100, 'r-', linewidth=1.5, label='Gene END', alpha=0.9)
    ax.axvline(0, color='black', linestyle='--', linewidth=0.8, alpha=0.5)
    ax.axhline(genome_at, color='gray', linestyle=':', linewidth=0.8, alpha=0.5, label=f'Genome AT ({genome_at:.1f}%)')
    ax.set_xlabel('Position relative to boundary (bp)')
    ax.set_ylabel('AT content (%)')
    ax.set_title('A. Meta-Gene AT Profile (Phase 9 Reproduced)')
    ax.legend(fontsize=7)
    ax.set_xlim(-500, 500)

    # Annotate the AT spike magnitudes
    # Find max in intergenic region (negative positions)
    ig_region = positions < 0
    if np.any(ig_region):
        start_spike = np.max(start_profile[ig_region]) * 100 - genome_at
        end_spike = np.max(end_profile[ig_region]) * 100 - genome_at
        ax.text(0.02, 0.95, f'Start spike: +{start_spike:.1f}pp\nEnd spike: +{end_spike:.1f}pp',
                transform=ax.transAxes, fontsize=7, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    # Panel B: Strand-aware boundary counts
    ax = axes[0, 1]
    labels = ['Promoter\nsignal', 'Terminator\nsignal']
    correct = [strand_data['promoter_signal'], strand_data['terminator_signal']]
    ax.bar(labels, correct, color=['royalblue', 'tomato'], edgecolor='black', alpha=0.8)
    ax.set_ylabel('Count (intergenic poly-A/T runs)')
    ax.set_title('B. Strand-Aware Regulatory Signal')

    total = strand_data.get('total_all_boundary', sum(correct))
    if total > 0:
        ax.text(0.02, 0.95, f'Total correct: {strand_data["total_correct"]}/{total}',
                transform=ax.transAxes, fontsize=8, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    # Panel C: Hairpin null model
    ax = axes[0, 2]
    obs = hairpin_data['observed_rate'] * 100
    null_mean = hairpin_data['null_mean'] * 100
    null_std = hairpin_data['null_std'] * 100
    excess = hairpin_data['excess'] * 100

    bar_labels = ['Observed', 'Null\n(Markov-1)']
    bar_vals = [obs, null_mean]
    bar_errs = [0, null_std]
    colors = ['#2ecc71', '#95a5a6']
    bars = ax.bar(bar_labels, bar_vals, yerr=bar_errs, color=colors,
                  edgecolor='black', alpha=0.8, capsize=5)
    ax.set_ylabel('Hairpin detection rate (%)')
    ax.set_title(f'C. Hairpin Null Model (stem>={hairpin_data["min_stem"]})')
    ax.text(0.02, 0.95,
            f'Excess: {excess:.1f}pp\nz={hairpin_data["z_score"]:.1f}, p={hairpin_data["p_value"]:.3f}',
            transform=ax.transAxes, fontsize=8, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    # Panel D: Two-population sizes
    ax = axes[1, 0]
    labels = ['Boundary\npoly-A', 'Boundary\npoly-T', 'Dispersed\npoly-A', 'Dispersed\npoly-T']
    vals = [twopop_data['n_boundary_A'], twopop_data['n_boundary_T'],
            twopop_data['n_dispersed_A'], twopop_data['n_dispersed_T']]
    colors = ['royalblue', 'tomato', 'lightblue', 'lightsalmon']
    ax.bar(labels, vals, color=colors, edgecolor='black', alpha=0.8)
    ax.set_ylabel('Count')
    ax.set_title('D. Two Populations: Boundary vs Dispersed')
    ax.text(0.02, 0.95, f'A/T ratio p={twopop_data["AT_ratio_p"]:.4f}',
            transform=ax.transAxes, fontsize=8, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='lightyellow', alpha=0.8))

    # Panel E: Plus strand detail
    ax = axes[1, 1]
    plus_labels = ['poly-A\nupstream\nof start\n(promoter)',
                   'poly-T\ndownstream\nof end\n(terminator)']
    plus_vals = [strand_data['pA_upstream_plus_start'],
                 strand_data['pT_downstream_plus_end']]
    ax.bar(plus_labels, plus_vals, color=['royalblue', 'tomato'],
           edgecolor='black', alpha=0.8)
    ax.set_ylabel('Count')
    ax.set_title('E. Plus Strand: Expected Signal Positions')

    # Panel F: Dispersed runs as latent regulatory substrate
    ax = axes[1, 2]
    if ig_poly_A is not None and ig_poly_T is not None and genes is not None:
        boundary_pos = []
        dispersed_pos = []

        starts = []
        ends_arr = []
        for g in genes:
            if g['strand'] == '+':
                starts.append(g['start_pos'])
                ends_arr.append(g['end_pos'])
            else:
                starts.append(g['end_pos'])
                ends_arr.append(g['start_pos'])
        starts = np.array(starts)
        ends_arr = np.array(ends_arr)

        all_ig_runs = ig_poly_A + ig_poly_T
        for r in all_ig_runs:
            pos = r['pos_1based']
            min_dist = min(np.min(np.abs(pos - starts)), np.min(np.abs(pos - ends_arr)))
            if min_dist <= 100:
                boundary_pos.append(pos)
            else:
                dispersed_pos.append(pos)

        bins = np.linspace(0, 2473870, 50)
        if boundary_pos:
            ax.hist(np.array(boundary_pos)/1e6, bins=bins/1e6, alpha=0.6,
                    color='green', label=f'Boundary ({len(boundary_pos)})', edgecolor='black', linewidth=0.5)
        if dispersed_pos:
            ax.hist(np.array(dispersed_pos)/1e6, bins=bins/1e6, alpha=0.4,
                    color='gray', label=f'Dispersed ({len(dispersed_pos)})', edgecolor='black', linewidth=0.5)
        ax.set_xlabel('Genome position (Mb)')
        ax.set_ylabel('Count')
        ax.set_title('F. Spatial Distribution: Active vs Latent')
        ax.legend(fontsize=7)
    else:
        ax.text(0.5, 0.5, 'Data not available', ha='center', va='center', transform=ax.transAxes)

    plt.tight_layout()
    outpath = FIGURES_DIR / 'phase1b_corrected.png'
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
    print("PANDORAVIRUS 2.0 -- PHASE 1b: CORRECTED ANALYSIS")
    print("=" * 70)

    conn = get_conn()
    cur = conn.cursor()

    # Load data
    print("\nLoading data...")
    seq = get_sequence(cur, genome_id=4)
    genes = get_genes(cur, genome_id=4)
    print(f"  Genome: {len(seq):,} bp, {len(genes)} genes")

    # Compute poly-A/T runs
    poly_A = find_poly_runs(seq, 'A', min_len=5)
    poly_T = find_poly_runs(seq, 'T', min_len=5)
    print(f"  poly-A runs (>=5): {len(poly_A)}")
    print(f"  poly-T runs (>=5): {len(poly_T)}")

    # ── FIX 1: Intergenic-only filter ──
    print("\n" + "=" * 70)
    print("FIX 1: Intergenic-Only Filter")
    print("=" * 70)

    ig_poly_A, cod_poly_A = filter_intergenic_runs(poly_A, genes)
    ig_poly_T, cod_poly_T = filter_intergenic_runs(poly_T, genes)

    print(f"  poly-A: {len(ig_poly_A)} intergenic, {len(cod_poly_A)} coding (removed)")
    print(f"  poly-T: {len(ig_poly_T)} intergenic, {len(cod_poly_T)} coding (removed)")
    print(f"  Total intergenic: {len(ig_poly_A) + len(ig_poly_T)} "
          f"(was {len(poly_A) + len(poly_T)}, removed {len(cod_poly_A) + len(cod_poly_T)})")

    # ── FIX 2: Strand-aware boundary test ──
    strand_data = strand_aware_boundary_analysis(seq, genes, ig_poly_A, ig_poly_T)

    # ── Meta-gene AT profile (Phase 9 reproduction) ──
    print("\n" + "=" * 70)
    print("META-GENE AT PROFILE (Phase 9 Reproduction)")
    print("=" * 70)
    positions, start_profile, end_profile = metagene_at_profile(seq, genes)

    genome_at = sum(1 for b in seq if b in 'AT') / len(seq)
    start_spike = (np.max(start_profile[positions < 0]) - genome_at) * 100
    end_spike = (np.max(end_profile[positions < 0]) - genome_at) * 100
    print(f"\n  Genome AT: {genome_at*100:.2f}%")
    print(f"  Gene START intergenic AT spike: +{start_spike:.1f} pp")
    print(f"  Gene END intergenic AT spike:   +{end_spike:.1f} pp")
    print(f"  (Phase 9 reported: +7.0pp start, +9.5pp end)")

    # ── FIX 3: Hairpin null model ──
    hairpin_data = hairpin_null_model(seq, genes, n_shuffles=20, min_stem=6)

    # ── Two-population test ──
    twopop_data = two_population_test(ig_poly_A, ig_poly_T, genes)

    # ── Visualization ──
    print("\n" + "=" * 70)
    print("GENERATING FIGURES")
    print("=" * 70)
    plot_phase1b_results(positions, start_profile, end_profile,
                         strand_data, hairpin_data, twopop_data,
                         ig_poly_A, ig_poly_T, genes)

    # ── Save results ──
    results = {
        'intergenic_filter': {
            'ig_poly_A': len(ig_poly_A),
            'ig_poly_T': len(ig_poly_T),
            'coding_removed_A': len(cod_poly_A),
            'coding_removed_T': len(cod_poly_T),
        },
        'strand_aware': strand_data,
        'metagene': {
            'start_spike_pp': start_spike,
            'end_spike_pp': end_spike,
        },
        'hairpin_null': hairpin_data,
        'two_population': twopop_data,
    }

    results_path = RESULTS_DIR / 'phase1b_results.json'
    with open(str(results_path), 'w') as f:
        json.dump(results, f, indent=2, default=str)
    print(f"  Results saved: {results_path}")

    # ── Summary ──
    elapsed = time.time() - t0
    print(f"\n" + "=" * 70)
    print(f"PHASE 1b CORRECTED SUMMARY")
    print(f"=" * 70)
    print(f"  Runtime: {elapsed:.1f}s")
    print(f"\n  FIX 1 (Intergenic filter):")
    print(f"    Removed {len(cod_poly_A) + len(cod_poly_T)} coding runs "
          f"({(len(cod_poly_A)+len(cod_poly_T))/(len(poly_A)+len(poly_T))*100:.1f}%)")
    print(f"\n  FIX 2 (Strand-aware boundary):")
    print(f"    Promoter signal: {strand_data['promoter_signal']}")
    print(f"    Terminator signal: {strand_data['terminator_signal']}")
    print(f"    Meta-gene AT spikes: +{start_spike:.1f}pp (start), +{end_spike:.1f}pp (end)")
    print(f"\n  FIX 3 (Hairpin null model):")
    print(f"    Observed: {hairpin_data['observed_rate']*100:.1f}%")
    print(f"    Null: {hairpin_data['null_mean']*100:.1f}% +/- {hairpin_data['null_std']*100:.1f}%")
    print(f"    Excess: {hairpin_data['excess']*100:.1f}pp (z={hairpin_data['z_score']:.1f})")
    print(f"\n  Two-population test:")
    print(f"    Boundary: {twopop_data['n_boundary_A'] + twopop_data['n_boundary_T']} runs")
    print(f"    Dispersed: {twopop_data['n_dispersed_A'] + twopop_data['n_dispersed_T']} runs")
    print(f"    A/T composition difference: p={twopop_data['AT_ratio_p']:.4f}")

    conn.close()


if __name__ == "__main__":
    main()
