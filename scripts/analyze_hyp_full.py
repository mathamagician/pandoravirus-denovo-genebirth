"""
Analyze Full Hypothetical BLAST Results + Update Paper
=======================================================
Run this AFTER blast_hyp_complete.py finishes.

What it does:
  1. Recomputes all hypothetical gene statistics using full N=902 dataset
  2. Computes new spatial OR with CI (should tighten from 2.32-42.83)
  3. Computes new annotated-vs-hypothetical Fisher OR with updated N
  4. Updates results/phase2_results.json with full-dataset figures
  5. Prints a summary of all numbers that need updating in draft4_GBE.md
  6. Saves results/hyp_full_analysis.json

Usage:
    python scripts/analyze_hyp_full.py
"""

import json
import math
import re
from pathlib import Path
from collections import Counter

import numpy as np
from scipy import stats

PROJECT_ROOT = Path(__file__).resolve().parent.parent
RESULTS_DIR  = PROJECT_ROOT / "results"

# ── Load data ──────────────────────────────────────────────────────────

def load_results():
    path = RESULTS_DIR / "blastp_hypothetical.json"
    with open(path) as f:
        return json.load(f)

def parse_pos(desc):
    m = re.search(r'pos:(\d+)-(\d+)', desc)
    return (int(m.group(1)), int(m.group(2))) if m else (None, None)

def haldane_or_ci(a, b, c, d, alpha=0.05):
    ah, bh, ch, dh = a+0.5, b+0.5, c+0.5, d+0.5
    ln_or = math.log(ah * dh / (bh * ch))
    se    = math.sqrt(1/ah + 1/bh + 1/ch + 1/dh)
    z     = 1.96
    return math.exp(ln_or - z*se), math.exp(ln_or), math.exp(ln_or + z*se)


def main():
    data = load_results()
    total_n = len(data)

    # Exclude genuine errors (no 'description' or network errors with no useful data)
    valid    = {k: v for k, v in data.items() if 'description' in v}
    errors   = {k: v for k, v in data.items() if 'error' in v}
    n_errors = len(errors)
    n_valid  = len(valid)

    # Count orphans vs hits
    n_hits    = sum(1 for v in valid.values() if v.get('n_hits_non_self', 0) > 0)
    n_orphans = n_valid - n_hits

    orphan_rate     = n_orphans / n_valid
    non_self_rate   = n_hits / n_valid

    print(f"\n{'='*60}")
    print(f"FULL HYPOTHETICAL BLAST ANALYSIS")
    print(f"{'='*60}")
    print(f"Total entries in JSON: {total_n}")
    print(f"Valid results:         {n_valid}")
    print(f"Errors (excluded):     {n_errors}")
    print(f"Non-self hits:         {n_hits} ({non_self_rate*100:.1f}%)")
    print(f"Orphans (no non-self): {n_orphans} ({orphan_rate*100:.1f}%)")

    # ── Spatial OR ───────────────────────────────────────────────────
    native_lo, native_hi = 1_300_000, 2_000_000

    hits_in_nz    = 0
    hits_out_nz   = 0
    nohits_in_nz  = 0
    nohits_out_nz = 0
    no_pos        = 0

    hit_organisms = []

    for v in valid.values():
        s, _ = parse_pos(v.get('description', ''))
        has_hit = v.get('n_hits_non_self', 0) > 0

        if s is None:
            no_pos += 1
            continue

        in_nz = native_lo <= s <= native_hi

        if in_nz:
            if has_hit: hits_in_nz += 1
            else:       nohits_in_nz += 1
        else:
            if has_hit: hits_out_nz += 1
            else:       nohits_out_nz += 1

        if has_hit:
            org = v.get('top_non_self_organism') or 'Unknown'
            hit_organisms.append(org)

    print(f"\n--- Spatial enrichment (native zone 1.3-2.0 Mb) ---")
    print(f"Entries without position: {no_pos}")
    print(f"2x2 table:")
    print(f"  hits_in_nz={hits_in_nz}, nohits_in_nz={nohits_in_nz}")
    print(f"  hits_out_nz={hits_out_nz}, nohits_out_nz={nohits_out_nz}")

    fisher_or, fisher_p = stats.fisher_exact(
        [[hits_in_nz, nohits_in_nz], [hits_out_nz, nohits_out_nz]]
    )
    ci_lo, hald_or, ci_hi = haldane_or_ci(hits_in_nz, nohits_in_nz, hits_out_nz, nohits_out_nz)

    print(f"Fisher OR = {fisher_or:.2f}, p = {fisher_p:.4f}")
    print(f"Haldane OR = {hald_or:.2f}, 95% CI: {ci_lo:.2f}–{ci_hi:.2f}")

    # ── Taxonomy of non-self hits ────────────────────────────────────
    print(f"\n--- Taxonomy of non-self hits (N={n_hits}) ---")
    org_counts = Counter(hit_organisms)
    for org, cnt in org_counts.most_common(15):
        print(f"  {cnt:3d}  {org}")

    # ── Updated annotated vs hypothetical comparison ─────────────────
    # Annotated numbers are fixed (from original analysis)
    ann_n    = 521
    ann_hits = 95
    hyp_n    = n_valid
    hyp_hits = n_hits

    ann_rate = ann_hits / ann_n
    hyp_rate = hyp_hits / hyp_n

    f_or, f_p = stats.fisher_exact([[ann_hits, ann_n - ann_hits], [hyp_hits, hyp_n - hyp_hits]])
    ci_lo2, hald_or2, ci_hi2 = haldane_or_ci(ann_hits, ann_n - ann_hits, hyp_hits, hyp_n - hyp_hits)

    print(f"\n--- Updated annotated vs hypothetical comparison ---")
    print(f"Annotated: {ann_hits}/{ann_n} = {ann_rate*100:.1f}% non-self hit rate")
    print(f"Hypothetical (full): {hyp_hits}/{hyp_n} = {hyp_rate*100:.1f}% non-self hit rate")
    print(f"Fisher OR = {f_or:.2f} (Haldane {hald_or2:.2f}, 95% CI: {ci_lo2:.2f}–{ci_hi2:.2f}), p = {f_p:.2e}")

    # ── Updated ORFan rate extrapolation ────────────────────────────
    # Now we have the full hypothetical count — no extrapolation needed
    total_genes    = 1430
    ann_orphans    = ann_n - ann_hits     # 426
    hyp_orphans    = n_orphans            # directly observed
    total_orphans  = ann_orphans + hyp_orphans
    # Use n_valid (successful queries) for denominator
    direct_orphan_rate = (ann_orphans + hyp_orphans) / (ann_n + hyp_n)

    print(f"\n--- Updated ORFan rate (direct, no extrapolation needed) ---")
    print(f"Annotated orphans: {ann_orphans}/{ann_n} = {ann_orphans/ann_n*100:.1f}%")
    print(f"Hypothetical orphans (full): {hyp_orphans}/{hyp_n} = {hyp_orphans/hyp_n*100:.1f}%")
    print(f"Combined direct rate: ({ann_orphans}+{hyp_orphans})/({ann_n}+{hyp_n}) = {direct_orphan_rate*100:.1f}%")
    print(f"Note: this is now a DIRECT observation, not an extrapolation")

    # ── Save results ─────────────────────────────────────────────────
    output = {
        "analysis_date": "2026-03-21",
        "total_hypothetical_queried": n_valid,
        "total_errors": n_errors,
        "non_self_hits": n_hits,
        "non_self_hit_rate": round(non_self_rate, 4),
        "orphans": n_orphans,
        "orphan_rate": round(orphan_rate, 4),
        "spatial_enrichment": {
            "hits_in_native_zone": hits_in_nz,
            "nohits_in_native_zone": nohits_in_nz,
            "hits_outside_native_zone": hits_out_nz,
            "nohits_outside_native_zone": nohits_out_nz,
            "fisher_or": round(float(fisher_or), 2),
            "fisher_p": round(float(fisher_p), 4),
            "haldane_or": round(hald_or, 2),
            "haldane_ci_lo": round(ci_lo, 2),
            "haldane_ci_hi": round(ci_hi, 2),
        },
        "annotated_vs_hypothetical": {
            "annotated_n": ann_n,
            "annotated_hits": ann_hits,
            "annotated_rate": round(ann_rate, 4),
            "hypothetical_n": hyp_n,
            "hypothetical_hits": hyp_hits,
            "hypothetical_rate": round(hyp_rate, 4),
            "fisher_or": round(float(f_or), 2),
            "haldane_or": round(hald_or2, 2),
            "haldane_ci_lo": round(ci_lo2, 2),
            "haldane_ci_hi": round(ci_hi2, 2),
            "fisher_p": float(f_p),
        },
        "direct_orphan_rate": {
            "annotated_orphans": int(ann_orphans),
            "hypothetical_orphans": int(hyp_orphans),
            "total_queried": int(ann_n + hyp_n),
            "combined_rate": round(direct_orphan_rate, 4),
            "note": "Direct observation from full dataset — no extrapolation needed"
        },
        "top_organisms": dict(org_counts.most_common(20)),
    }

    out_path = RESULTS_DIR / "hyp_full_analysis.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved to {out_path}")

    # ── Paper update checklist ───────────────────────────────────────
    print(f"\n{'='*60}")
    print("NUMBERS TO UPDATE IN PAPER (draft4 → draft5)")
    print(f"{'='*60}")
    print(f"Abstract:")
    print(f"  hypothetical orphan rate: 97.0% → {orphan_rate*100:.1f}%  (now direct, not sampled)")
    print(f"  overall orphan rate: 91.4% (extrapolated) → {direct_orphan_rate*100:.1f}% (direct)")
    print(f"Results (homology section):")
    print(f"  '299 successfully queried hypothetical' → '{hyp_n} successfully queried'")
    print(f"  '3.0%' non-self hit rate → '{non_self_rate*100:.1f}%'")
    print(f"  OR = 7.19 (annotated vs hyp) → OR = {f_or:.2f} (Haldane {hald_or2:.2f}, CI {ci_lo2:.2f}–{ci_hi2:.2f})")
    print(f"Spatial OR (hypothetical):")
    print(f"  Old: OR=11.7, CI 2.32–42.83  →  New: OR={fisher_or:.2f}, Haldane {hald_or:.2f}, CI {ci_lo:.2f}–{ci_hi:.2f}")
    print(f"Limitations:")
    print(f"  Remove: 'The hypothetical gene subsample (n=299) represents approximately one-third...'")
    print(f"  Replace with: 'All 902 hypothetical genes were queried; N={hyp_n} returned valid results.'")


if __name__ == "__main__":
    main()
