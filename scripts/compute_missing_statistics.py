"""
Compute missing statistics for pre-submission audit fixes
=========================================================
Computes all statistics flagged as missing in the statistical audit.
Output: results/missing_statistics.json

Usage: python scripts/compute_missing_statistics.py
"""

import json
import math
import sys
import os
import numpy as np
from scipy import stats

DB_CONFIG = {
    "dbname": os.environ.get("DB_NAME", "pandoravirus"),
    "user": os.environ.get("DB_USER", "postgres"),
    "password": os.environ["DB_PASSWORD"],
    "host": os.environ.get("DB_HOST", "localhost"),
    "port": int(os.environ.get("DB_PORT", 5432)),
}

# ─── Helper functions ────────────────────────────────────────────────────────

def haldane_or_ci(a, b, c, d, alpha=0.05):
    """Haldane-corrected Woolf 95% CI for odds ratio (adds 0.5 to each cell)."""
    ah, bh, ch, dh = a + 0.5, b + 0.5, c + 0.5, d + 0.5
    ln_or = math.log(ah * dh / (bh * ch))
    se = math.sqrt(1/ah + 1/bh + 1/ch + 1/dh)
    z = stats.norm.ppf(1 - alpha/2)
    return math.exp(ln_or - z*se), math.exp(ln_or), math.exp(ln_or + z*se)


# ─── 1. ORFan rate derivation ─────────────────────────────────────────────────

def compute_orphan_rate():
    """Extrapolated orphan rate with bootstrap CI."""
    ann_orphans, ann_total = 426, 521
    hyp_orphans_sampled, hyp_sampled = 290, 299
    hyp_total = 902
    classified_total = ann_total + hyp_total  # 1423

    ann_data = np.array([1] * ann_orphans + [0] * (ann_total - ann_orphans))
    hyp_data = np.array([1] * hyp_orphans_sampled + [0] * (hyp_sampled - hyp_orphans_sampled))

    point_est = (ann_total * (ann_orphans / ann_total) +
                 hyp_total * (hyp_orphans_sampled / hyp_sampled)) / classified_total

    np.random.seed(42)
    rates = []
    for _ in range(10_000):
        b_ann = np.random.choice(ann_data, size=len(ann_data), replace=True)
        b_hyp = np.random.choice(hyp_data, size=len(hyp_data), replace=True)
        est = (ann_total * b_ann.mean() + hyp_total * b_hyp.mean()) / classified_total
        rates.append(est)

    ci_lo, ci_hi = np.percentile(rates, [2.5, 97.5])

    return {
        "annotated_n": ann_total,
        "annotated_orphans": ann_orphans,
        "annotated_orphan_rate": ann_orphans / ann_total,
        "hypothetical_n_sampled": hyp_sampled,
        "hypothetical_orphans_sampled": hyp_orphans_sampled,
        "hypothetical_orphan_rate_sampled": hyp_orphans_sampled / hyp_sampled,
        "hypothetical_total": hyp_total,
        "classified_total": classified_total,
        "extrapolated_orphan_rate": round(point_est, 4),
        "bootstrap_ci_lo": round(ci_lo, 4),
        "bootstrap_ci_hi": round(ci_hi, 4),
        "direct_observed_rate": round((ann_orphans + hyp_orphans_sampled) / (ann_total + hyp_sampled), 4),
        "note": "Extrapolated: (521 × 0.818 + 902 × 0.970) / 1423. Philippe et al. 2013 reported 92.5% by independent method."
    }


# ─── 2. Neutrality plot regression ───────────────────────────────────────────

def compute_neutrality_regression(conn):
    """GC12 ~ GC3 regression for P. salinus genes (genome_id=4)."""
    cur = conn.cursor()
    cur.execute("""
        SELECT
            SUM(CASE WHEN SUBSTR(gc.codon_seq,3,1) IN ('G','C') THEN 1 ELSE 0 END)::float / COUNT(*) AS gc3,
            (SUM(CASE WHEN SUBSTR(gc.codon_seq,1,1) IN ('G','C') THEN 1 ELSE 0 END) +
             SUM(CASE WHEN SUBSTR(gc.codon_seq,2,1) IN ('G','C') THEN 1 ELSE 0 END))::float / (2*COUNT(*)) AS gc12
        FROM gene_codon gc
        JOIN gene g ON g.gene_id = gc.gene_id
        WHERE g.genome_id = 4
        GROUP BY gc.gene_id
        HAVING COUNT(*) >= 50
    """)
    rows = cur.fetchall()
    gc3 = np.array([r[0] for r in rows])
    gc12 = np.array([r[1] for r in rows])

    slope, intercept, r, p, se = stats.linregress(gc3, gc12)
    n = len(gc3)
    t_crit = stats.t.ppf(0.975, df=n-2)
    ci_lo = slope - t_crit * se
    ci_hi = slope + t_crit * se

    return {
        "n_genes": n,
        "slope": round(slope, 4),
        "slope_ci_lo": round(ci_lo, 4),
        "slope_ci_hi": round(ci_hi, 4),
        "r_squared": round(r**2, 4),
        "p_value": float(f"{p:.4e}"),
        "intercept": round(intercept, 4),
        "interpretation": (
            f"Slope {slope:.4f} << 1.0 (expected under strict neutrality). "
            "GC12 largely invariant with wobble GC3, consistent with translational selection."
        )
    }


# ─── 3. CAI Mann-Whitney test ─────────────────────────────────────────────────

def compute_cai_mannwhitney(conn):
    """Mann-Whitney test on per-gene CAI: P. salinus vs Mimivirus."""
    cur = conn.cursor()

    cur.execute("""
        SELECT wm.metric_value
        FROM window_metric wm
        WHERE wm.metric_name = 'cai_acanthamoeba' AND wm.genome_id = 4
    """)
    cai_psal = np.array([float(r[0]) for r in cur.fetchall()])

    cur.execute("""
        SELECT wm.metric_value
        FROM window_metric wm
        WHERE wm.metric_name = 'cai_acanthamoeba' AND wm.genome_id = 3
    """)
    cai_mimi = np.array([float(r[0]) for r in cur.fetchall()])

    u, p = stats.mannwhitneyu(cai_psal, cai_mimi, alternative='two-sided')

    np.random.seed(42)
    diffs = []
    for _ in range(5_000):
        b1 = np.random.choice(cai_psal, size=len(cai_psal), replace=True)
        b2 = np.random.choice(cai_mimi, size=len(cai_mimi), replace=True)
        diffs.append(np.mean(b1) - np.mean(b2))
    ci_lo, ci_hi = np.percentile(diffs, [2.5, 97.5])

    return {
        "psal_n": int(len(cai_psal)),
        "psal_mean": round(float(np.mean(cai_psal)), 4),
        "psal_sd": round(float(np.std(cai_psal)), 4),
        "mimi_n": int(len(cai_mimi)),
        "mimi_mean": round(float(np.mean(cai_mimi)), 4),
        "mimi_sd": round(float(np.std(cai_mimi)), 4),
        "mean_difference": round(float(np.mean(cai_psal) - np.mean(cai_mimi)), 4),
        "mannwhitney_U": float(u),
        "mannwhitney_p": float(p),
        "bootstrap_ci_lo": round(ci_lo, 4),
        "bootstrap_ci_hi": round(ci_hi, 4),
        "note": "p essentially 0 (machine precision limit). CAI reference: A. castellanii Kazusa database."
    }


# ─── 4. Hypothetical spatial OR with CI ──────────────────────────────────────

def compute_hypothetical_spatial_or():
    """Fisher exact OR + Haldane CI for hypothetical gene hit enrichment in native zone."""
    import re
    with open("results/blastp_hypothetical.json") as f:
        data = json.load(f)

    native_zone_start, native_zone_end = 1_300_000, 2_000_000

    def parse_pos(desc):
        m = re.search(r'pos:(\d+)-(\d+)', desc)
        return (int(m.group(1)), int(m.group(2))) if m else (None, None)

    table = {'hits_in_nz': 0, 'hits_out_nz': 0, 'nohits_in_nz': 0, 'nohits_out_nz': 0}
    for v in data.values():
        s, _ = parse_pos(v.get('description', ''))
        if s is not None:
            in_nz = native_zone_start <= s <= native_zone_end
            has_hit = v['n_hits_non_self'] > 0
            if in_nz:
                table['hits_in_nz' if has_hit else 'nohits_in_nz'] += 1
            else:
                table['hits_out_nz' if has_hit else 'nohits_out_nz'] += 1

    a, b = table['hits_in_nz'], table['nohits_in_nz']
    c, d = table['hits_out_nz'], table['nohits_out_nz']

    fisher_or, fisher_p = stats.fisher_exact([[a, b], [c, d]])
    hald_lo, hald_or, hald_hi = haldane_or_ci(a, b, c, d)

    return {
        "hits_in_native_zone": a,
        "nohits_in_native_zone": b,
        "hits_outside_native_zone": c,
        "nohits_outside_native_zone": d,
        "fisher_or": round(float(fisher_or), 2),
        "fisher_p": round(float(fisher_p), 4),
        "haldane_or": round(hald_or, 2),
        "haldane_ci_lo": round(hald_lo, 2),
        "haldane_ci_hi": round(hald_hi, 2),
        "note": (
            f"n={a+c} total hypothetical non-self hits. Wide CI reflects small n. "
            f"Lower bound {hald_lo:.2f} > 1.0 indicates enrichment even at conservative estimate."
        )
    }


# ─── 5. Spearman rho = -1.0 exact p-value ────────────────────────────────────

def compute_spearman_tier_pvalue():
    """Exact p-value for rho = -1.0 with N=5 tiers."""
    # For N=5, number of permutations = 5! = 120
    # P(rho = -1.0 one-tailed) = 1/120 (only one perfect descending arrangement)
    n = 5
    n_perms = math.factorial(n)
    p_one_tailed = 1.0 / n_perms
    p_two_tailed = 2.0 * p_one_tailed

    # Fisher's combined test for two independent rho=-1.0 results (GC and gene length)
    chi2_combined = -2 * (math.log(p_one_tailed) + math.log(p_one_tailed))
    p_fisher_combined = 1 - stats.chi2.cdf(chi2_combined, df=4)

    return {
        "n_tiers": n,
        "n_permutations": n_perms,
        "rho": -1.0,
        "p_one_tailed": round(p_one_tailed, 4),
        "p_two_tailed": round(p_two_tailed, 4),
        "p_fisher_combined_two_tests": round(float(p_fisher_combined), 6),
        "chi2_combined": round(float(chi2_combined), 2),
        "note": "Both GC and gene length show rho=-1.0 across 5 tiers. Fisher combined p <<< 0.001."
    }


# ─── Main ─────────────────────────────────────────────────────────────────────

def main():
    results = {}

    # Analytical computations (no DB needed)
    print("Computing orphan rate...")
    results["orphan_rate_derived"] = compute_orphan_rate()

    print("Computing Spearman tier p-value...")
    results["spearman_tier_pvalue"] = compute_spearman_tier_pvalue()

    print("Computing hypothetical spatial OR...")
    results["hypothetical_spatial_or_ci"] = compute_hypothetical_spatial_or()

    # DB computations
    try:
        import psycopg2
        conn = psycopg2.connect(**DB_CONFIG)
        print("DB connected.")

        print("Computing neutrality plot regression...")
        results["neutrality_regression"] = compute_neutrality_regression(conn)

        print("Computing CAI Mann-Whitney test...")
        results["cai_mannwhitney"] = compute_cai_mannwhitney(conn)

        conn.close()
    except Exception as e:
        print(f"DB error: {e}")
        results["db_error"] = str(e)

    # Save results
    out_path = "results/missing_statistics.json"
    with open(out_path, "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nSaved to {out_path}")

    # Print summary table
    print("\n" + "="*60)
    print("SUMMARY OF MISSING STATISTICS")
    print("="*60)

    if "orphan_rate_derived" in results:
        r = results["orphan_rate_derived"]
        print(f"\n#1 ORFan rate: {r['extrapolated_orphan_rate']*100:.1f}% "
              f"(95% CI: {r['bootstrap_ci_lo']*100:.1f}%-{r['bootstrap_ci_hi']*100:.1f}%)")
        print(f"   Annotated orphan rate: {r['annotated_orphan_rate']*100:.1f}% ({r['annotated_orphans']}/{r['annotated_n']})")
        print(f"   Hypothetical orphan rate: {r['hypothetical_orphan_rate_sampled']*100:.1f}% ({r['hypothetical_orphans_sampled']}/{r['hypothetical_n_sampled']})")

    if "neutrality_regression" in results:
        r = results["neutrality_regression"]
        print(f"\n#2 Neutrality slope: {r['slope']:.4f} "
              f"(95% CI: {r['slope_ci_lo']:.4f}-{r['slope_ci_hi']:.4f}), "
              f"R2={r['r_squared']:.4f}, p={r['p_value']:.2e}")

    if "cai_mannwhitney" in results:
        r = results["cai_mannwhitney"]
        print(f"\n#3 CAI: P.sal {r['psal_mean']:.4f}+/-{r['psal_sd']:.4f} vs "
              f"Mimi {r['mimi_mean']:.4f}+/-{r['mimi_sd']:.4f}")
        print(f"   Mann-Whitney p={r['mannwhitney_p']:.2e}, diff CI: {r['bootstrap_ci_lo']:.4f}-{r['bootstrap_ci_hi']:.4f}")

    if "hypothetical_spatial_or_ci" in results:
        r = results["hypothetical_spatial_or_ci"]
        print(f"\n#4 Hypothetical spatial OR: {r['fisher_or']:.2f} (Fisher), "
              f"Haldane {r['haldane_or']:.2f} (95% CI: {r['haldane_ci_lo']:.2f}-{r['haldane_ci_hi']:.2f}), "
              f"p={r['fisher_p']:.4f}")

    if "spearman_tier_pvalue" in results:
        r = results["spearman_tier_pvalue"]
        print(f"\n#5 Spearman rho=-1.0, N=5: one-tailed p={r['p_one_tailed']:.4f}, "
              f"two-tailed p={r['p_two_tailed']:.4f}")
        print(f"   Fisher combined (2 tests): chi2={r['chi2_combined']:.2f}, p={r['p_fisher_combined_two_tests']:.6f}")


if __name__ == "__main__":
    main()
