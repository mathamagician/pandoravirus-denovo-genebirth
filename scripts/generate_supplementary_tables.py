"""
Generate all supplementary tables (S1-S9) for the combined paper.
Pulls data from PostgreSQL database and results JSON files.
Outputs CSV files to supplementary/ directory.
"""

import json
import csv
import os
import psycopg2
import psycopg2.extras

# Config
DB_PARAMS = dict(host='localhost', port=5432, dbname='pandoravirus', user='postgres', password='pandora2026')
RESULTS_DIR = os.path.join(os.path.dirname(__file__), '..', 'results')
OUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'supplementary')
os.makedirs(OUT_DIR, exist_ok=True)

def get_conn():
    return psycopg2.connect(**DB_PARAMS)


def table_s1_codon_usage():
    """Table S1: Codon usage comparison between annotated and ORFan genes."""
    conn = get_conn()
    cur = conn.cursor()

    # Get per-gene codon metrics, split by annotation status
    # window_metric uses gene_id as window_start for per-gene metrics
    cur.execute("""
        SELECT
            g.product,
            g.gc_content,
            wm_gc3.metric_value AS gc3,
            wm_gc12.metric_value AS gc12,
            wm_enc.metric_value AS enc,
            wm_cai.metric_value AS cai
        FROM gene g
        LEFT JOIN window_metric wm_gc3 ON wm_gc3.genome_id = g.genome_id
            AND wm_gc3.metric_name = 'gc3' AND wm_gc3.window_start = g.gene_id
        LEFT JOIN window_metric wm_gc12 ON wm_gc12.genome_id = g.genome_id
            AND wm_gc12.metric_name = 'gc12' AND wm_gc12.window_start = g.gene_id
        LEFT JOIN window_metric wm_enc ON wm_enc.genome_id = g.genome_id
            AND wm_enc.metric_name = 'enc' AND wm_enc.window_start = g.gene_id
        LEFT JOIN window_metric wm_cai ON wm_cai.genome_id = g.genome_id
            AND wm_cai.metric_name = 'cai_acanthamoeba' AND wm_cai.window_start = g.gene_id
        WHERE g.genome_id = 4
    """)

    rows = cur.fetchall()
    conn.close()

    # Classify: annotated vs hypothetical/ORFan
    annotated = []
    orfan = []
    for product, gc, gc3, gc12, enc, cai in rows:
        entry = {'gc': float(gc) if gc else None,
                 'gc3': float(gc3) if gc3 else None,
                 'gc12': float(gc12) if gc12 else None,
                 'enc': float(enc) if enc else None,
                 'cai': float(cai) if cai else None}
        if product and 'hypothetical' not in product.lower():
            annotated.append(entry)
        else:
            orfan.append(entry)

    def mean_of(lst, key):
        vals = [x[key] for x in lst if x[key] is not None]
        return sum(vals) / len(vals) if vals else None

    def std_of(lst, key):
        vals = [x[key] for x in lst if x[key] is not None]
        if len(vals) < 2:
            return None
        m = sum(vals) / len(vals)
        return (sum((v - m)**2 for v in vals) / (len(vals) - 1)) ** 0.5

    metrics = ['gc', 'gc3', 'gc12', 'enc', 'cai']
    metric_labels = ['GC content', 'GC3 (wobble position)', 'GC12 (constrained positions)',
                     'ENC (effective number of codons)', 'CAI (codon adaptation index)']

    with open(os.path.join(OUT_DIR, 'table_s1_codon_usage.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['Metric', 'Annotated_n', 'Annotated_mean', 'Annotated_sd',
                     'ORFan_n', 'ORFan_mean', 'ORFan_sd', 'Difference'])
        for metric, label in zip(metrics, metric_labels):
            a_vals = [x[metric] for x in annotated if x[metric] is not None]
            o_vals = [x[metric] for x in orfan if x[metric] is not None]
            a_mean = mean_of(annotated, metric)
            o_mean = mean_of(orfan, metric)
            diff = (a_mean - o_mean) if (a_mean and o_mean) else None
            w.writerow([label, len(a_vals),
                        f'{a_mean:.4f}' if a_mean else 'N/A',
                        f'{std_of(annotated, metric):.4f}' if std_of(annotated, metric) else 'N/A',
                        len(o_vals),
                        f'{o_mean:.4f}' if o_mean else 'N/A',
                        f'{std_of(orfan, metric):.4f}' if std_of(orfan, metric) else 'N/A',
                        f'{diff:+.4f}' if diff else 'N/A'])

    print(f"  S1: {len(annotated)} annotated, {len(orfan)} ORFan genes")


def table_s2_structural_families():
    """Table S2: ESM-2 structural family assignments (from Phase 8 gene_classification)."""
    conn = get_conn()
    cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)

    # Get Phase 8 structural classifications
    cur.execute("""
        SELECT gc.gene_id, gc.classification, gc.confidence, gc.notes,
               g.product, g.gc_content, g.gene_length, g.start_pos, g.end_pos, g.strand
        FROM gene_classification gc
        JOIN gene g ON gc.gene_id = g.gene_id
        WHERE gc.source = 'structure_phase8'
        ORDER BY gc.classification, g.start_pos
    """)
    rows = cur.fetchall()
    conn.close()

    with open(os.path.join(OUT_DIR, 'table_s2_structural_families.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['Gene_ID', 'Structural_family', 'Confidence', 'Product',
                     'GC_content', 'Gene_length_bp', 'Start_pos', 'End_pos', 'Strand', 'Notes'])
        for row in rows:
            w.writerow([row['gene_id'], row['classification'],
                        f"{row['confidence']:.3f}" if row['confidence'] else 'N/A',
                        row['product'] or 'hypothetical protein',
                        f"{row['gc_content']:.4f}" if row['gc_content'] else 'N/A',
                        row['gene_length'], row['start_pos'], row['end_pos'],
                        row['strand'], row['notes'] or ''])

    print(f"  S2: {len(rows)} proteins with structural family assignments")


def table_s3_gene_families():
    """Table S3: Gene family detection results (from Phase 7 gene_classification)."""
    conn = get_conn()
    cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)

    cur.execute("""
        SELECT gc.gene_id, gc.classification, gc.confidence, gc.notes,
               g.product, g.gc_content, g.gene_length, g.start_pos
        FROM gene_classification gc
        JOIN gene g ON gc.gene_id = g.gene_id
        WHERE gc.source = 'gene_family_phase7'
        ORDER BY gc.classification, g.start_pos
    """)
    rows = cur.fetchall()
    conn.close()

    with open(os.path.join(OUT_DIR, 'table_s3_gene_families.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['Gene_ID', 'Family', 'Similarity_score', 'Product',
                     'GC_content', 'Gene_length_bp', 'Genomic_position_bp'])
        for row in rows:
            w.writerow([row['gene_id'], row['classification'],
                        f"{row['confidence']:.3f}" if row['confidence'] else 'N/A',
                        row['product'] or 'hypothetical protein',
                        f"{row['gc_content']:.4f}" if row['gc_content'] else 'N/A',
                        row['gene_length'], row['start_pos']])

    print(f"  S3: {len(rows)} genes in families")


def table_s4_product_hit_rates():
    """Table S4: Product-level non-self BLASTp hit rates."""
    with open(os.path.join(RESULTS_DIR, 'phase2_results.json')) as f:
        data = json.load(f)

    products = data.get('products', {})

    with open(os.path.join(OUT_DIR, 'table_s4_product_hit_rates.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['Product', 'Total_genes', 'Non_self_hits', 'Hit_rate_pct', 'Top_hit_categories'])
        for product, info in sorted(products.items(), key=lambda x: -x[1].get('hit_rate', 0)):
            w.writerow([product, info.get('total', ''), info.get('non_self', ''),
                        f"{info.get('hit_rate', 0):.1f}" if 'hit_rate' in info else 'N/A',
                        info.get('top_categories', '')])

    print(f"  S4: {len(products)} product categories")


def table_s5_evalue_sensitivity():
    """Table S5: E-value sensitivity analysis."""
    with open(os.path.join(RESULTS_DIR, 'phase2_results.json')) as f:
        data = json.load(f)

    sensitivity = data.get('sensitivity', {})

    with open(os.path.join(OUT_DIR, 'table_s5_evalue_sensitivity.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['E_value_threshold', 'Non_self_hits', 'Hit_rate_pct', 'Orphan_rate_pct'])
        for threshold, info in sorted(sensitivity.items(), key=lambda x: float(x[0])):
            hits = info.get('non_self_hits', info.get('hits', ''))
            hit_rate = info.get('hit_rate', '')
            orphan_rate = info.get('orphan_rate', '')
            w.writerow([threshold, hits,
                        f"{hit_rate:.1f}" if isinstance(hit_rate, (int, float)) else hit_rate,
                        f"{orphan_rate:.1f}" if isinstance(orphan_rate, (int, float)) else orphan_rate])

    print(f"  S5: {len(sensitivity)} E-value thresholds")


def table_s6_all_nonself_hits():
    """Table S6: Complete taxonomy of all 95 non-self BLAST hits."""
    with open(os.path.join(RESULTS_DIR, 'phase2_results.json')) as f:
        data = json.load(f)

    protein_tax = data.get('protein_taxonomy', {})

    # Get gene positions from database
    conn = get_conn()
    cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
    cur.execute("SELECT protein_id, product, start_pos, end_pos, strand FROM gene WHERE genome_id = 4")
    gene_info = {row['protein_id']: row for row in cur.fetchall()}
    conn.close()

    with open(os.path.join(OUT_DIR, 'table_s6_nonself_hits_annotated.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['Protein_ID', 'Product', 'Genomic_position_bp', 'Strand',
                     'Taxonomy_category', 'Top_hit_organism', 'Pct_identity',
                     'E_value', 'Query_coverage', 'In_native_zone'])

        count = 0
        for pid, info in sorted(protein_tax.items()):
            if info['category'] == 'No non-self hit':
                continue
            gi = gene_info.get(pid, {})
            start = gi.get('start_pos', '')
            in_zone = 'Yes' if (isinstance(start, int) and 1300000 <= start <= 2000000) else 'No'
            w.writerow([pid, gi.get('product', ''), start, gi.get('strand', ''),
                        info['category'], info['organism'],
                        f"{info['pct_identity']:.1f}" if info['pct_identity'] else '',
                        f"{info['evalue']:.2e}" if info['evalue'] else '',
                        f"{info.get('query_coverage', '')}" if info.get('query_coverage') else '',
                        in_zone])
            count += 1

    print(f"  S6: {count} non-self hits for annotated genes")


def table_s7_hypothetical_hits():
    """Table S7: All 9 hypothetical gene non-self hits."""
    with open(os.path.join(RESULTS_DIR, 'blastp_hypothetical.json')) as f:
        hyp_data = json.load(f)  # dict keyed by protein_id

    # Get gene positions
    conn = get_conn()
    cur = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
    cur.execute("SELECT protein_id, start_pos, end_pos, strand FROM gene WHERE genome_id = 4")
    gene_info = {row['protein_id']: row for row in cur.fetchall()}
    conn.close()

    with open(os.path.join(OUT_DIR, 'table_s7_nonself_hits_hypothetical.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        w.writerow(['Protein_ID', 'Genomic_position_bp', 'Strand',
                     'Taxonomy_category', 'Top_hit_organism', 'Pct_identity',
                     'E_value', 'In_native_zone'])

        count = 0
        for pid, entry in hyp_data.items():
            if entry.get('n_hits_non_self', 0) == 0:
                continue

            hits = entry.get('hits', [])

            # Find first non-self hit
            non_self_hit = None
            for hit in hits:
                org = hit.get('organism', '')
                if org and 'pandora' not in org.lower():
                    non_self_hit = hit
                    break

            if non_self_hit is None:
                continue

            gi = gene_info.get(pid, {})
            start = gi.get('start_pos', '')
            in_zone = 'Yes' if (isinstance(start, int) and 1300000 <= start <= 2000000) else 'No'

            org_name = non_self_hit.get('organism', 'Unknown')
            identity = non_self_hit.get('identity', non_self_hit.get('pct_identity', ''))
            evalue = non_self_hit.get('evalue', non_self_hit.get('e_value', ''))

            # Classify taxonomy
            org_lower = org_name.lower() if org_name else ''
            if any(x in org_lower for x in ['cedrat', 'molli', 'bodo', 'medusa', 'mimivirus', 'marseillevirus', 'faust', 'pitho', 'kaumoeba']):
                cat = 'Giant virus'
            elif any(x in org_lower for x in ['acanthamoeba', 'balamuthia', 'amoeb']):
                cat = 'Amoebozoa'
            elif any(x in org_lower for x in ['micromonas', 'bathycoccus', 'ostreococcus', 'ectocarpus', 'scytosiphon']):
                cat = 'Marine algae'
            else:
                cat = 'Bacteria'

            w.writerow([pid, start, gi.get('strand', ''),
                        cat, org_name,
                        f"{identity:.1f}" if isinstance(identity, (int, float)) else identity,
                        f"{evalue:.2e}" if isinstance(evalue, float) else evalue,
                        in_zone])
            count += 1

    print(f"  S7: {count} non-self hits for hypothetical genes")


def table_s8_crossgenome_cv():
    """Table S8: Cross-genome compositional homogeneity (full metrics)."""
    conn = get_conn()
    cur = conn.cursor()

    genomes = {1: 'PhiX174', 2: 'Lambda', 3: 'Mimivirus', 4: 'P. salinus'}
    metrics = ['gc_content', 'gc3', 'gc12', 'enc', 'cai_acanthamoeba']
    metric_labels = ['GC', 'GC3', 'GC12', 'ENC', 'CAI']

    results = {}
    for gid, gname in genomes.items():
        results[gname] = {}
        for metric in metrics:
            if metric == 'gc_content':
                cur.execute("""
                    SELECT AVG(gc_content), STDDEV(gc_content), COUNT(*)
                    FROM gene WHERE genome_id = %s AND gc_content IS NOT NULL
                """, (gid,))
            else:
                cur.execute("""
                    SELECT AVG(metric_value), STDDEV(metric_value), COUNT(*)
                    FROM window_metric
                    WHERE genome_id = %s AND metric_name = %s AND metric_value IS NOT NULL
                """, (gid, metric))

            avg, std, n = cur.fetchone()
            if avg and std and avg != 0:
                cv = (std / avg) * 100
            else:
                cv = None
            results[gname][metric] = {'mean': avg, 'sd': std, 'n': n, 'cv': cv}

    conn.close()

    with open(os.path.join(OUT_DIR, 'table_s8_crossgenome_cv.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        header = ['Genome']
        for label in metric_labels:
            header.extend([f'{label}_mean', f'{label}_sd', f'{label}_CV_pct', f'{label}_n'])
        w.writerow(header)

        for gname in ['PhiX174', 'Lambda', 'Mimivirus', 'P. salinus']:
            row = [gname]
            for metric in metrics:
                d = results[gname][metric]
                row.append(f"{d['mean']:.4f}" if d['mean'] else 'N/A')
                row.append(f"{d['sd']:.4f}" if d['sd'] else 'N/A')
                row.append(f"{d['cv']:.1f}" if d['cv'] else 'N/A')
                row.append(d['n'] or 0)
            w.writerow(row)

    print(f"  S8: Cross-genome CV for {len(genomes)} genomes x {len(metrics)} metrics")


def table_s9_hypothesis_matrix():
    """Table S9: Hypothesis testing matrix across all 12 modules."""

    matrix = [
        ['Module', 'Finding', 'HGT', 'Gene_duplication', 'Mobile_elements', 'De_novo_birth'],
        ['2: Genome architecture', 'GC gap 10.7pp, lowest coding density (67.5%)',
         'Neutral', 'Neutral', 'Neutral', 'Consistent (large intergenic space = substrate)'],
        ['3: Null models', 'Poly-A/T 15x enriched; AGCT 1500-fold depleted; 88.6% effect-size noise',
         'Neutral', 'Neutral', 'Neutral', 'Consistent (regulatory signals genuine, not artifacts)'],
        ['4: Codon usage', 'ORFan GC3 diff +0.015, CAI diff +0.008; upstream cosine 0.9824',
         'Disfavored (no foreign signatures)', 'Neutral', 'Neutral', 'Consistent (compositional homogeneity)'],
        ['5: Gene prediction', '347 Prodigal-only genes; true IG drops to 23.8%',
         'Neutral', 'Neutral', 'Neutral', 'Consistent (unannotated proto-genes exist)'],
        ['6: Evolutionary classification', 'GC CV 6.1% (lowest); 3 GMM clusters; no foreign outliers',
         'Disfavored (no compositional heterogeneity)', 'Neutral', 'Neutral', 'Consistent (uniform composition)'],
        ['7: Gene families', '96.8% singletons; 14 domain-repeat families; no tandem clustering; null d=0.06',
         'Neutral', 'Disfavored (no whole-gene duplication)', 'Disfavored (no expansion signal)', 'Consistent (singletons expected)'],
        ['8: Protein structure', '9 structural families; d=1.45 vs null; 41% GC confound; proto-gene NN p=0.90',
         'Neutral', 'Neutral', 'Neutral', 'Consistent (structural continuity across tiers)'],
        ['7b-c: Proto-gene continuum', 'T1-T5 monotonic gradient; T5 codon p=0.0005; no spatial clustering',
         'Disfavored (no sharp boundaries)', 'Neutral', 'Neutral', 'Supported (gradual maturation)'],
        ['9: Regulatory signals', '+10.3/+9.8pp AT spikes; 82.8% intergenic; T4 p=0.14, T5 p=0.25',
         'Neutral', 'Neutral', 'Neutral', 'Supported (distributed regulatory infrastructure)'],
        ['10: Transcriptomic validation', '74% strand-aware correct (p<10^-6); two-pop AT p=0.105',
         'Neutral', 'Neutral', 'Neutral', 'Supported (regulatory signals validated)'],
        ['11-12: BLASTp homology', '81.8% orphan; OR=7.19 (p=7.79e-12); phylogenetic distance gradient',
         'Disfavored (eukaryote-dominated, not multi-kingdom)', 'Neutral', 'Neutral', 'Supported (age-dependent homology gradient)'],
    ]

    with open(os.path.join(OUT_DIR, 'table_s9_hypothesis_matrix.csv'), 'w', newline='') as f:
        w = csv.writer(f)
        for row in matrix:
            w.writerow(row)

    print(f"  S9: Hypothesis matrix with {len(matrix)-1} modules x 4 hypotheses")


if __name__ == '__main__':
    print("Generating supplementary tables...")
    print()

    try:
        table_s1_codon_usage()
    except Exception as e:
        print(f"  S1 ERROR: {e}")

    try:
        table_s2_structural_families()
    except Exception as e:
        print(f"  S2 ERROR: {e}")

    try:
        table_s3_gene_families()
    except Exception as e:
        print(f"  S3 ERROR: {e}")

    try:
        table_s4_product_hit_rates()
    except Exception as e:
        print(f"  S4 ERROR: {e}")

    try:
        table_s5_evalue_sensitivity()
    except Exception as e:
        print(f"  S5 ERROR: {e}")

    try:
        table_s6_all_nonself_hits()
    except Exception as e:
        print(f"  S6 ERROR: {e}")

    try:
        table_s7_hypothetical_hits()
    except Exception as e:
        print(f"  S7 ERROR: {e}")

    try:
        table_s8_crossgenome_cv()
    except Exception as e:
        print(f"  S8 ERROR: {e}")

    try:
        table_s9_hypothesis_matrix()
    except Exception as e:
        print(f"  S9 ERROR: {e}")

    print()
    print("Done. Files written to supplementary/ directory.")
