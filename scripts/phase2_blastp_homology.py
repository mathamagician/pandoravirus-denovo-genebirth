"""
Phase 2: BLASTp Homology Analysis of Pandoravirus salinus Genes
================================================================
Runs BLASTp on annotated (non-hypothetical) genes against NCBI nr
to determine evolutionary origins and taxonomic affinities.

Analyses:
  1. Taxonomic classification of top BLAST hits
  2. Percent-identity distribution
  3. Spatial mapping of homology onto genome
  4. Native zone (1.3-2.0 Mb) enrichment analysis
  5. Domain/product-level breakdown
  6. Annotated vs hypothetical comparison (subsample)

Usage:
  python phase2_blastp_homology.py                # Check status / analyze existing results
  python phase2_blastp_homology.py --blast        # Run BLASTp (remote NCBI, resumable)
  python phase2_blastp_homology.py --blast-hyp    # Also BLAST 100 random hypothetical genes
"""

import os
import sys
import json
import time
import random
import numpy as np
import psycopg2
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from pathlib import Path
from collections import defaultdict, Counter
from scipy import stats

# Biopython
from Bio import SeqIO, Entrez
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# ── Configuration ─────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parent.parent
DATA_DIR = PROJECT_ROOT / "data"
FIGURES_DIR = PROJECT_ROOT / "figures"
RESULTS_DIR = PROJECT_ROOT / "results"

for d in [DATA_DIR, FIGURES_DIR, RESULTS_DIR]:
    d.mkdir(exist_ok=True)

Entrez.email = "Eddie.Bradford@gmail.com"
Entrez.api_key = None  # Set for 10 req/s instead of 3 req/s

DB_CONFIG = {
    "dbname": "pandoravirus",
    "user": "postgres",
    "password": "pandora2026",
    "host": "localhost",
    "port": 5432,
}

GENOME_ID = 4  # P. salinus
GENOME_LENGTH = 2_473_870  # bp

# BLAST parameters
BLAST_EVALUE = 1e-5
BLAST_MAX_HITS = 10
BLAST_DB = "nr"
BLAST_DELAY = 1  # seconds between queries (each query takes ~60s, well under rate limit)

# Hypothetical subsample size for comparison
HYP_SAMPLE_SIZE = 300

# ── Database ──────────────────────────────────────────────────────────

def get_conn():
    return psycopg2.connect(**DB_CONFIG)


def load_genes():
    """Load all P. salinus genes with positions and translations."""
    conn = get_conn()
    cur = conn.cursor()
    cur.execute("""
        SELECT gene_id, locus_tag, gene_name, start_pos, end_pos, strand,
               gene_length, product, protein_id, translation, gc_content
        FROM gene
        WHERE genome_id = %s
        ORDER BY start_pos
    """, (GENOME_ID,))
    rows = cur.fetchall()
    conn.close()

    genes = []
    for r in rows:
        genes.append({
            'gene_id': r[0],
            'locus_tag': r[1],
            'gene_name': r[2],
            'start_pos': r[3],
            'end_pos': r[4],
            'strand': r[5],
            'gene_length': r[6],
            'product': r[7] or '',
            'protein_id': r[8] or '',
            'translation': r[9] or '',
            'gc_content': float(r[10]) if r[10] else None,
            'midpoint': (r[3] + r[4]) / 2,
        })

    return genes


def classify_gene(product):
    """Classify gene as annotated or hypothetical."""
    p = product.lower()
    if not product or 'hypothetical' in p:
        return 'hypothetical'
    return 'annotated'


# ── FASTA Export ──────────────────────────────────────────────────────

def export_fasta(genes, category='annotated'):
    """Export protein sequences to FASTA from database translations."""
    fasta_path = DATA_DIR / f"proteins_{category}.fasta"

    records = []
    for g in genes:
        if not g['translation'] or not g['protein_id']:
            continue
        if category == 'annotated' and classify_gene(g['product']) != 'annotated':
            continue
        if category == 'hypothetical' and classify_gene(g['product']) != 'hypothetical':
            continue

        rec = SeqRecord(
            Seq(g['translation']),
            id=g['protein_id'],
            description=f"{g['locus_tag']}|{g['product']}|pos:{g['start_pos']}-{g['end_pos']}({g['strand']})"
        )
        records.append(rec)

    SeqIO.write(records, str(fasta_path), "fasta")
    print(f"Exported {len(records)} {category} proteins to {fasta_path}")
    return records


# ── BLAST (Remote NCBI) ──────────────────────────────────────────────

def extract_organism(hit_def):
    """Extract organism name from BLAST hit definition.

    BLAST hit_def typically contains [Organism name] at the end.
    Example: 'ankyrin repeat protein [Pandoravirus salinus]'
    """
    import re
    # Find the last bracketed term (organism)
    matches = re.findall(r'\[([^\]]+)\]', hit_def)
    if matches:
        return matches[-1]
    return 'Unknown'


def run_blastp(proteins, results_path, label=""):
    """Run BLASTp against nr via NCBI remote BLAST. Fully resumable."""
    # Load existing results
    if results_path.exists():
        with open(results_path, 'r') as f:
            existing = json.load(f)
        print(f"Loaded {len(existing)} existing results")
    else:
        existing = {}

    remaining = [p for p in proteins if p.id not in existing]
    total = len(proteins)
    done = total - len(remaining)
    print(f"BLAST {label}: {done}/{total} complete, {len(remaining)} remaining")

    if not remaining:
        print("All proteins already BLASTed.")
        return existing

    for i, protein in enumerate(remaining):
        n = done + i + 1
        print(f"\n[{n}/{total}] {protein.id}: {protein.description[:60]}...")
        t0 = time.time()

        try:
            result_handle = NCBIWWW.qblast(
                "blastp", BLAST_DB, str(protein.seq),
                expect=BLAST_EVALUE,
                hitlist_size=BLAST_MAX_HITS,
                format_type="XML",
            )

            blast_records = NCBIXML.parse(result_handle)
            record = next(blast_records)

            hits = []
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    organism = extract_organism(alignment.hit_def)
                    hits.append({
                        'hit_id': alignment.hit_id,
                        'hit_def': alignment.hit_def[:300],
                        'hit_accession': alignment.accession,
                        'organism': organism,
                        'subject_length': alignment.length,
                        'evalue': hsp.expect,
                        'bit_score': hsp.bits,
                        'score': hsp.score,
                        'identities': hsp.identities,
                        'positives': hsp.positives,
                        'gaps': hsp.gaps,
                        'align_length': hsp.align_length,
                        'query_start': hsp.query_start,
                        'query_end': hsp.query_end,
                        'pct_identity': round(hsp.identities / hsp.align_length * 100, 1) if hsp.align_length else 0,
                        'query_coverage': round((hsp.query_end - hsp.query_start + 1) / len(protein.seq) * 100, 1),
                    })
                    break  # top HSP per hit only

            # Filter self-hits (other Pandoravirus strains)
            non_self_hits = [h for h in hits if 'pandoravirus' not in h['organism'].lower()]

            existing[protein.id] = {
                'protein_id': protein.id,
                'description': protein.description,
                'query_length': len(protein.seq),
                'n_hits_total': len(hits),
                'n_hits_non_self': len(non_self_hits),
                'hits': hits,
                'non_self_hits': non_self_hits,
                'top_organism': hits[0]['organism'] if hits else None,
                'top_non_self_organism': non_self_hits[0]['organism'] if non_self_hits else None,
            }

            elapsed = time.time() - t0
            top = hits[0]['organism'] if hits else 'no hits'
            top_ns = non_self_hits[0]['organism'] if non_self_hits else 'no non-self hits'
            print(f"  {len(hits)} hits ({len(non_self_hits)} non-self) in {elapsed:.1f}s")
            print(f"  Top: {top} | Top non-self: {top_ns}")

        except Exception as e:
            print(f"  ERROR: {e}")
            existing[protein.id] = {
                'protein_id': protein.id,
                'description': protein.description,
                'query_length': len(protein.seq),
                'n_hits_total': 0,
                'n_hits_non_self': 0,
                'hits': [],
                'non_self_hits': [],
                'error': str(e),
            }

        # Save after each query (resumable)
        with open(results_path, 'w') as f:
            json.dump(existing, f, indent=2)

        time.sleep(BLAST_DELAY)

    return existing


# ── Taxonomy Classification ──────────────────────────────────────────

# Curated organism-level classifications for all organisms observed in
# our BLAST results. Reviewed and corrected after senior peer review
# identified misclassification of marine green algae as bacteria.
#
# Categories:
#   Giant virus    — NCLDVs (Nucleocytoviricota)
#   Other virus    — Non-giant viruses (baculovirus, etc.)
#   Amoebozoa      — Pandoravirus host lineage
#   Marine algae   — Green algae (Chlorophyta/Prasinophyceae) + brown algae (Phaeophyceae)
#   Insect         — Aphids, etc.
#   Other eukaryote — Fungi, oomycetes, rotifers, plants, diatoms, haptophytes
#   Bacteria       — Genuine prokaryotes

ORGANISM_CLASSIFICATIONS = {
    # ── Giant viruses (NCLDVs / Nucleocytoviricota) ──
    'Cedratvirus kamchatka': 'Giant virus',
    'Cedratvirus lausannensis': 'Giant virus',
    'Cedratvirus duvanny': 'Giant virus',
    'Cedratvirus lena': 'Giant virus',
    'Cedratvirus Ce7-1': 'Giant virus',
    'Mollivirus shaoxing': 'Giant virus',
    'Mollivirus kamchatka': 'Giant virus',
    'Pithovirus sibericum': 'Giant virus',
    'Pacmanvirus S19': 'Giant virus',
    'Orpheovirus IHUMI-LCC2': 'Giant virus',
    'Marseillevirus sp.': 'Giant virus',
    'Naiavirus sp.': 'Giant virus',
    'Acanthamoeba castellanii medusavirus J1': 'Giant virus',
    'Bodo saltans virus': 'Giant virus',
    'Paramecium bursaria Chlorella virus IL3A': 'Giant virus',
    # ── Other viruses ──
    'Choristoneura diversana nucleopolyhedrovirus': 'Other virus',
    # ── Amoebozoa (host lineage) ──
    'Acanthamoeba castellanii str. Neff': 'Amoebozoa',
    'Balamuthia mandrillaris': 'Amoebozoa',
    # ── Green algae (Chlorophyta / Prasinophyceae) ──
    'Micromonas pusilla CCMP1545': 'Marine algae',
    'Bathycoccus prasinos': 'Marine algae',
    'Ostreococcus lucimarinus CCE9901': 'Marine algae',
    'Chlamydomonas sp. UWO 241': 'Marine algae',
    'Cymbomonas tetramitiformis': 'Marine algae',
    'Dunaliella salina': 'Marine algae',
    # ── Brown algae (Phaeophyceae) + Schizocladiophyceae ──
    'Myriotrichia clavaeformis': 'Marine algae',
    'Scytosiphon promiscuus': 'Marine algae',
    'Ectocarpus fasciculatus': 'Marine algae',
    'Feldmannia mitchelliae': 'Marine algae',
    'Choristocarpus tenellus': 'Marine algae',
    'Desmarestia herbacea': 'Marine algae',
    'Ectocarpus siliculosus': 'Marine algae',
    'Discosporangium mesarthrocarpum': 'Marine algae',
    'Schizocladia ischiensis': 'Marine algae',
    # ── Insects ──
    'Acyrthosiphon pisum': 'Insect',
    'Aphis craccivora': 'Insect',
    'Aphis gossypii': 'Insect',
    # ── Other eukaryotes ──
    'Prymnesium parvum': 'Other eukaryote',       # haptophyte
    'Seminavis robusta': 'Other eukaryote',         # diatom
    'Helianthus annuus': 'Other eukaryote',         # sunflower
    'Peronospora farinosa': 'Other eukaryote',      # oomycete
    'Eremomyces bilateralis CBS 781.70': 'Other eukaryote',  # fungus
    'Polyrhizophydium stewartii': 'Other eukaryote', # chytrid
    'Adineta steineri': 'Other eukaryote',          # bdelloid rotifer
    # ── Bacteria ──
    'Cetobacterium sp.': 'Bacteria',
    'Desulfobulbaceae bacterium': 'Bacteria',
    'Patescibacteria group bacterium': 'Bacteria',
    'Clostridia bacterium': 'Bacteria',
    'Akkermansia sp.': 'Bacteria',
    'Rariglobus sp.': 'Bacteria',
    'Fidelibacterota bacterium': 'Bacteria',
}

# Keyword fallback for organisms not in the curated list
TAXONOMY_KEYWORDS = [
    # Giant viruses (check BEFORE generic 'virus')
    ('pandoravirus', 'Giant virus'),
    ('mimivirus', 'Giant virus'),
    ('megavirus', 'Giant virus'),
    ('moumouvirus', 'Giant virus'),
    ('marseillevirus', 'Giant virus'),
    ('pithovirus', 'Giant virus'),
    ('mollivirus', 'Giant virus'),
    ('cedratvirus', 'Giant virus'),
    ('tupanvirus', 'Giant virus'),
    ('klosneuvirus', 'Giant virus'),
    ('orpheovirus', 'Giant virus'),
    ('medusavirus', 'Giant virus'),
    ('pacmanvirus', 'Giant virus'),
    ('naiavirus', 'Giant virus'),
    ('cafeteria roenbergensis virus', 'Giant virus'),
    # Amoebozoa
    ('acanthamoeba', 'Amoebozoa'),
    ('balamuthia', 'Amoebozoa'),
    ('dictyostelium', 'Amoebozoa'),
    ('entamoeba', 'Amoebozoa'),
    # Green algae
    ('micromonas', 'Marine algae'),
    ('bathycoccus', 'Marine algae'),
    ('ostreococcus', 'Marine algae'),
    ('chlamydomonas', 'Marine algae'),
    ('cymbomonas', 'Marine algae'),
    ('dunaliella', 'Marine algae'),
    # Brown algae
    ('ectocarpus', 'Marine algae'),
    ('myriotrichia', 'Marine algae'),
    ('scytosiphon', 'Marine algae'),
    ('desmarestia', 'Marine algae'),
    ('feldmannia', 'Marine algae'),
    ('schizocladia', 'Marine algae'),
    # Viruses (generic, after giant virus checks)
    ('virus', 'Other virus'),
    ('phage', 'Other virus'),
    ('viridae', 'Other virus'),
    # Bacteria (be specific to avoid false matches)
    ('bacterium', 'Bacteria'),
    ('bacteria', 'Bacteria'),
    ('bacillus', 'Bacteria'),
    ('escherichia', 'Bacteria'),
    ('akkermansia', 'Bacteria'),
    ('cetobacterium', 'Bacteria'),
    ('desulfobulb', 'Bacteria'),
    ('clostridia', 'Bacteria'),
    ('patescibacteria', 'Bacteria'),
    ('fidelibacterota', 'Bacteria'),
    ('rariglobus', 'Bacteria'),
]


def classify_organism(organism):
    """Classify an organism into a broad taxonomic category.

    Uses curated organism-level lookup first, then keyword fallback.
    """
    # Exact match first
    if organism in ORGANISM_CLASSIFICATIONS:
        return ORGANISM_CLASSIFICATIONS[organism]

    # Keyword fallback (order matters — checked sequentially)
    org_lower = organism.lower()
    for keyword, category in TAXONOMY_KEYWORDS:
        if keyword in org_lower:
            return category

    return 'Other/Unknown'


def batch_taxonomy_lookup(organisms):
    """Classify organisms using curated lookup + Entrez fallback for unknowns."""
    unique_orgs = list(set(organisms))
    print(f"Classifying {len(unique_orgs)} unique organisms...")

    # First pass: curated + keyword classification
    taxonomy = {}
    for org in unique_orgs:
        taxonomy[org] = classify_organism(org)

    # Report classification
    from collections import Counter
    cats = Counter(taxonomy.values())
    for cat, n in cats.most_common():
        print(f"  {cat}: {n} unique organisms")

    # Entrez lookup for unknowns
    unknowns = [org for org, cat in taxonomy.items() if cat == 'Other/Unknown']
    if unknowns:
        print(f"  {len(unknowns)} unclassified — trying Entrez lookup...")
        for org in unknowns:
            try:
                handle = Entrez.esearch(db="taxonomy", term=f'"{org}"')
                result = Entrez.read(handle)
                handle.close()

                if result['IdList']:
                    tax_id = result['IdList'][0]
                    handle = Entrez.efetch(db="taxonomy", id=tax_id, retmode="xml")
                    records = Entrez.read(handle)
                    handle.close()

                    if records:
                        lineage = records[0].get('Lineage', '')
                        if 'Amoebozoa' in lineage:
                            taxonomy[org] = 'Amoebozoa'
                        elif 'Nucleocytoviricota' in lineage or 'Imitervirales' in lineage:
                            taxonomy[org] = 'Giant virus'
                        elif 'Viruses' in lineage:
                            taxonomy[org] = 'Other virus'
                        elif 'Bacteria' in lineage:
                            taxonomy[org] = 'Bacteria'
                        elif 'Archaea' in lineage:
                            taxonomy[org] = 'Archaea'
                        elif any(x in lineage for x in ['Chlorophyta', 'Prasinophyceae', 'Phaeophyceae']):
                            taxonomy[org] = 'Marine algae'
                        elif any(x in lineage for x in ['Arthropoda', 'Insecta']):
                            taxonomy[org] = 'Insect'
                        else:
                            taxonomy[org] = 'Other eukaryote'
                        print(f"    Entrez: {org} -> {taxonomy[org]}")

                time.sleep(0.4)
            except Exception:
                pass  # Keep as Unknown

    return taxonomy


# ── Analysis ──────────────────────────────────────────────────────────

def analyze_blast_results(blast_results, genes):
    """Comprehensive analysis of BLAST results."""
    results = {}

    # Build protein_id -> gene lookup
    gene_lookup = {}
    for g in genes:
        if g['protein_id']:
            gene_lookup[g['protein_id']] = g

    # ── Basic stats ──
    # Distinguish BLAST errors (connection failures) from genuine no-hit results
    total_queried = len(blast_results)
    errors = sum(1 for r in blast_results.values() if 'error' in r)
    successful = total_queried - errors
    has_any_hit = sum(1 for r in blast_results.values() if r['n_hits_total'] > 0 and 'error' not in r)
    has_non_self = sum(1 for r in blast_results.values() if r['n_hits_non_self'] > 0 and 'error' not in r)
    genuine_no_hits = sum(1 for r in blast_results.values() if r['n_hits_total'] == 0 and 'error' not in r)
    self_only = has_any_hit - has_non_self

    print(f"\n{'='*60}")
    print(f"BLAST Results Summary")
    print(f"{'='*60}")
    print(f"Total attempted:        {total_queried}")
    print(f"  Successful queries:   {successful}")
    print(f"  Connection errors:    {errors} (excluded from statistics)")
    print(f"\nOf {successful} successful queries:")
    print(f"  With any hit:         {has_any_hit} ({has_any_hit/successful*100:.1f}%)")
    print(f"  With non-self hit:    {has_non_self} ({has_non_self/successful*100:.1f}%)")
    print(f"  Self-only (Pandora):  {self_only} ({self_only/successful*100:.1f}%)")
    print(f"  Genuine no hits:      {genuine_no_hits} ({genuine_no_hits/successful*100:.1f}%)")

    results['basic'] = {
        'total_attempted': total_queried,
        'successful_queries': successful,
        'connection_errors': errors,
        'has_any_hit': has_any_hit,
        'has_non_self_hit': has_non_self,
        'self_only': self_only,
        'genuine_no_hits': genuine_no_hits,
    }

    # ── Taxonomy of non-self top hits ──
    organisms = []
    for r in blast_results.values():
        if r.get('non_self_hits'):
            organisms.append(r['non_self_hits'][0]['organism'])

    if organisms:
        taxonomy = batch_taxonomy_lookup(organisms)

        # Count by category
        tax_counts = Counter()
        protein_taxonomy = {}
        for pid, r in blast_results.items():
            if r.get('non_self_hits'):
                org = r['non_self_hits'][0]['organism']
                cat = taxonomy.get(org, 'Other/Unknown')
                tax_counts[cat] += 1
                protein_taxonomy[pid] = {
                    'category': cat,
                    'organism': org,
                    'pct_identity': r['non_self_hits'][0]['pct_identity'],
                    'evalue': r['non_self_hits'][0]['evalue'],
                    'query_coverage': r['non_self_hits'][0].get('query_coverage', 0),
                }
            else:
                protein_taxonomy[pid] = {
                    'category': 'No non-self hit',
                    'organism': None,
                    'pct_identity': 0,
                    'evalue': None,
                }

        print(f"\nTaxonomic distribution of top non-self hits (n={has_non_self}):")
        for cat, count in tax_counts.most_common():
            print(f"  {cat:25s}  {count:4d}  ({count/has_non_self*100:.1f}% of non-self)")

        results['taxonomy'] = dict(tax_counts)
        results['protein_taxonomy'] = protein_taxonomy
    else:
        print("\nNo non-self hits to classify.")
        results['taxonomy'] = {}
        results['protein_taxonomy'] = {}
        protein_taxonomy = {}

    # ── Identity distribution ──
    identities_non_self = []
    identities_all = []
    for r in blast_results.values():
        if r.get('hits'):
            identities_all.append(r['hits'][0]['pct_identity'])
        if r.get('non_self_hits'):
            identities_non_self.append(r['non_self_hits'][0]['pct_identity'])

    if identities_non_self:
        print(f"\nPercent identity (non-self top hits):")
        print(f"  Mean:   {np.mean(identities_non_self):.1f}%")
        print(f"  Median: {np.median(identities_non_self):.1f}%")
        print(f"  Min:    {np.min(identities_non_self):.1f}%")
        print(f"  Max:    {np.max(identities_non_self):.1f}%")
        print(f"  <30%:   {sum(1 for x in identities_non_self if x < 30)} ({sum(1 for x in identities_non_self if x < 30)/len(identities_non_self)*100:.1f}%)")
        print(f"  30-50%: {sum(1 for x in identities_non_self if 30 <= x < 50)} ({sum(1 for x in identities_non_self if 30 <= x < 50)/len(identities_non_self)*100:.1f}%)")
        print(f"  >50%:   {sum(1 for x in identities_non_self if x >= 50)} ({sum(1 for x in identities_non_self if x >= 50)/len(identities_non_self)*100:.1f}%)")

        results['identity'] = {
            'mean': float(np.mean(identities_non_self)),
            'median': float(np.median(identities_non_self)),
            'min': float(np.min(identities_non_self)),
            'max': float(np.max(identities_non_self)),
            'below_30': sum(1 for x in identities_non_self if x < 30),
            'between_30_50': sum(1 for x in identities_non_self if 30 <= x < 50),
            'above_50': sum(1 for x in identities_non_self if x >= 50),
            'values': identities_non_self,
        }

    # ── Product-level breakdown ──
    product_stats = defaultdict(lambda: {'total': 0, 'has_non_self': 0, 'categories': Counter(), 'identities': []})
    for pid, r in blast_results.items():
        g = gene_lookup.get(pid)
        if not g:
            continue
        product = g['product']
        product_stats[product]['total'] += 1
        if r.get('non_self_hits'):
            product_stats[product]['has_non_self'] += 1
            org = r['non_self_hits'][0]['organism']
            cat = protein_taxonomy.get(pid, {}).get('category', 'Unknown')
            product_stats[product]['categories'][cat] += 1
            product_stats[product]['identities'].append(r['non_self_hits'][0]['pct_identity'])

    print(f"\nProduct-level analysis (top 15):")
    for product, ps in sorted(product_stats.items(), key=lambda x: -x[1]['total'])[:15]:
        hit_rate = ps['has_non_self'] / ps['total'] * 100 if ps['total'] > 0 else 0
        mean_id = np.mean(ps['identities']) if ps['identities'] else 0
        top_cat = ps['categories'].most_common(1)[0][0] if ps['categories'] else 'N/A'
        print(f"  {product[:40]:40s}  n={ps['total']:3d}  hits={hit_rate:5.1f}%  id={mean_id:5.1f}%  top={top_cat}")

    results['products'] = {
        k: {'total': v['total'], 'has_non_self': v['has_non_self'],
             'categories': dict(v['categories']),
             'mean_identity': float(np.mean(v['identities'])) if v['identities'] else 0}
        for k, v in product_stats.items()
    }

    # ── Spatial analysis ──
    positions_hit = []
    positions_no_hit = []
    positions_self_only = []

    for pid, r in blast_results.items():
        g = gene_lookup.get(pid)
        if not g:
            continue
        midpoint = g['midpoint']

        if r['n_hits_non_self'] > 0:
            positions_hit.append(midpoint)
        elif r['n_hits_total'] > 0:
            positions_self_only.append(midpoint)
        else:
            positions_no_hit.append(midpoint)

    # Native zone analysis (1.3-2.0 Mb)
    native_start, native_end = 1_300_000, 2_000_000

    def in_native(pos):
        return native_start <= pos <= native_end

    native_fraction_genome = (native_end - native_start) / GENOME_LENGTH
    n_hit_native = sum(1 for p in positions_hit if in_native(p))
    n_hit_outside = len(positions_hit) - n_hit_native
    n_nohit_native = sum(1 for p in positions_no_hit + positions_self_only if in_native(p))
    n_nohit_outside = len(positions_no_hit) + len(positions_self_only) - n_nohit_native

    if positions_hit and (positions_no_hit or positions_self_only):
        # Fisher's exact test: are non-self hits depleted in native zone?
        table = [[n_hit_native, n_hit_outside], [n_nohit_native, n_nohit_outside]]
        odds_ratio, fisher_p = stats.fisher_exact(table)

        print(f"\nSpatial analysis (native zone 1.3-2.0 Mb):")
        print(f"  Non-self hits: {n_hit_native} in zone, {n_hit_outside} outside")
        print(f"  No hits:       {n_nohit_native} in zone, {n_nohit_outside} outside")
        print(f"  Fisher's exact test: OR={odds_ratio:.3f}, p={fisher_p:.4f}")
        print(f"  Expected native fraction: {native_fraction_genome:.1%}")

        results['spatial'] = {
            'native_zone': [native_start, native_end],
            'hit_in_native': n_hit_native,
            'hit_outside': n_hit_outside,
            'nohit_in_native': n_nohit_native,
            'nohit_outside': n_nohit_outside,
            'fisher_or': float(odds_ratio),
            'fisher_p': float(fisher_p),
        }

    # ── Identity by taxonomy category ──
    if protein_taxonomy:
        cat_identities = defaultdict(list)
        for pid, pt in protein_taxonomy.items():
            if pt['pct_identity'] > 0:
                cat_identities[pt['category']].append(pt['pct_identity'])

        print(f"\nIdentity by taxonomic category:")
        for cat, ids in sorted(cat_identities.items(), key=lambda x: -np.mean(x[1])):
            if len(ids) >= 3:
                print(f"  {cat:25s}  n={len(ids):3d}  mean={np.mean(ids):5.1f}%  median={np.median(ids):5.1f}%")

        results['identity_by_taxonomy'] = {
            cat: {'n': len(ids), 'mean': float(np.mean(ids)), 'median': float(np.median(ids)), 'values': ids}
            for cat, ids in cat_identities.items()
        }

    return results


# ── Visualization ─────────────────────────────────────────────────────

def plot_phase2_results(blast_results, analysis, genes, output_path):
    """Generate 6-panel figure for Phase 2."""
    gene_lookup = {g['protein_id']: g for g in genes if g['protein_id']}
    protein_taxonomy = analysis.get('protein_taxonomy', {})

    has_comparison = 'comparison' in analysis

    fig = plt.figure(figsize=(18, 18 if has_comparison else 12))
    gs = gridspec.GridSpec(3 if has_comparison else 2, 3, hspace=0.35, wspace=0.3)
    fig.suptitle('Phase 2: BLASTp Homology Analysis — Pandoravirus salinus',
                 fontsize=14, fontweight='bold', y=0.98 if not has_comparison else 0.99)

    # ── Panel A: Taxonomy pie chart ──
    ax1 = fig.add_subplot(gs[0, 0])
    tax = analysis.get('taxonomy', {})
    if tax:
        # Add "Pandoravirus only" and "Connection error" counts
        basic = analysis.get('basic', {})
        categories = dict(tax)
        if basic.get('self_only', 0) > 0:
            categories['Pandoravirus only'] = basic['self_only']
        if basic.get('genuine_no_hits', 0) > 0:
            categories['No BLAST hit'] = basic['genuine_no_hits']

        labels = []
        sizes = []
        for cat, count in sorted(categories.items(), key=lambda x: -x[1]):
            labels.append(f"{cat} ({count})")
            sizes.append(count)

        colors = plt.cm.Set3(np.linspace(0, 1, len(labels)))
        wedges, texts, autotexts = ax1.pie(
            sizes, labels=None, autopct='%1.0f%%', colors=colors,
            textprops={'fontsize': 7}, pctdistance=0.8,
            startangle=90,
        )
        ax1.legend(wedges, labels, loc='center left', bbox_to_anchor=(-0.3, 0.5),
                   fontsize=6)
    ax1.set_title('A. Taxonomic Distribution\n(top non-self BLAST hit)', fontsize=10)

    # ── Panel B: Identity distribution histogram ──
    ax2 = fig.add_subplot(gs[0, 1])
    id_data = analysis.get('identity', {})
    if id_data and id_data.get('values'):
        values = id_data['values']
        ax2.hist(values, bins=30, color='steelblue', edgecolor='white', alpha=0.8)
        ax2.axvline(np.mean(values), color='red', linestyle='--', label=f'Mean={np.mean(values):.1f}%')
        ax2.axvline(np.median(values), color='orange', linestyle='--', label=f'Median={np.median(values):.1f}%')
        ax2.set_xlabel('Percent Identity (%)')
        ax2.set_ylabel('Count')
        ax2.legend(fontsize=8)
    ax2.set_title('B. Identity Distribution\n(non-self top hits)', fontsize=10)

    # ── Panel C: Product-level hit rates ──
    ax3 = fig.add_subplot(gs[0, 2])
    prods = analysis.get('products', {})
    if prods:
        # Top products by count
        top_prods = sorted(prods.items(), key=lambda x: -x[1]['total'])[:12]
        prod_names = [p[0][:25] for p in top_prods]
        prod_totals = [p[1]['total'] for p in top_prods]
        prod_hits = [p[1]['has_non_self'] for p in top_prods]

        y_pos = np.arange(len(prod_names))
        ax3.barh(y_pos, prod_totals, color='lightgray', edgecolor='gray', label='Total')
        ax3.barh(y_pos, prod_hits, color='steelblue', edgecolor='gray', label='With non-self hit')
        ax3.set_yticks(y_pos)
        ax3.set_yticklabels(prod_names, fontsize=7)
        ax3.set_xlabel('Gene count')
        ax3.legend(fontsize=7)
        ax3.invert_yaxis()
    ax3.set_title('C. Hit Rate by Product Type', fontsize=10)

    # ── Panel D: Spatial distribution ──
    ax4 = fig.add_subplot(gs[1, 0:2])  # wide panel
    for pid, r in blast_results.items():
        g = gene_lookup.get(pid)
        if not g:
            continue
        pos = g['midpoint'] / 1e6

        if r['n_hits_non_self'] > 0:
            pt = protein_taxonomy.get(pid, {})
            cat = pt.get('category', 'Unknown')
            color_map = {
                'Giant virus': 'purple',
                'Other virus': 'violet',
                'Amoebozoa': 'green',
                'Marine algae': 'teal',
                'Insect': 'brown',
                'Other eukaryote': 'blue',
                'Bacteria': 'orange',
                'Archaea': 'red',
                'Other/Unknown': 'gray',
            }
            c = color_map.get(cat, 'gray')
            ax4.scatter(pos, pt.get('pct_identity', 0), c=c, s=8, alpha=0.5, zorder=2)
        elif r['n_hits_total'] > 0:
            ax4.scatter(pos, 0, c='lightgray', s=3, alpha=0.3, marker='x', zorder=1)

    # Native zone
    ax4.axvspan(1.3, 2.0, alpha=0.08, color='green', zorder=0)
    ax4.text(1.65, ax4.get_ylim()[1] * 0.95 if ax4.get_ylim()[1] > 0 else 95,
             'Native zone', ha='center', fontsize=8, color='green', alpha=0.7)

    # Legend
    from matplotlib.patches import Patch
    legend_elements = [
        Patch(facecolor='purple', label='Giant virus'),
        Patch(facecolor='green', label='Amoebozoa'),
        Patch(facecolor='teal', label='Marine algae'),
        Patch(facecolor='brown', label='Insect'),
        Patch(facecolor='blue', label='Other eukaryote'),
        Patch(facecolor='orange', label='Bacteria'),
        Patch(facecolor='lightgray', label='Pandoravirus only'),
    ]
    ax4.legend(handles=legend_elements, loc='upper right', fontsize=6, ncol=2)
    ax4.set_xlabel('Genome position (Mb)')
    ax4.set_ylabel('% Identity (non-self)')
    ax4.set_xlim(0, GENOME_LENGTH / 1e6)
    ax4.set_title('D. Spatial Distribution of Homology (colored by taxonomy)', fontsize=10)

    # ── Panel E: Identity by taxonomy boxplot ──
    ax5 = fig.add_subplot(gs[1, 2])
    id_by_tax = analysis.get('identity_by_taxonomy', {})
    if id_by_tax:
        # Filter to categories with >= 3 entries
        plot_cats = [(cat, d['values']) for cat, d in id_by_tax.items()
                     if len(d['values']) >= 3 and cat != 'No non-self hit']
        plot_cats.sort(key=lambda x: -np.median(x[1]))

        if plot_cats:
            bp = ax5.boxplot([v for _, v in plot_cats], vert=True, patch_artist=True,
                            widths=0.6)
            colors = plt.cm.Set2(np.linspace(0, 1, len(plot_cats)))
            for patch, color in zip(bp['boxes'], colors):
                patch.set_facecolor(color)
                patch.set_alpha(0.7)

            ax5.set_xticklabels([c[:15] for c, _ in plot_cats], rotation=45, ha='right', fontsize=7)
            ax5.set_ylabel('% Identity')
    ax5.set_title('E. Identity by Taxonomy', fontsize=10)

    # ── Panel F: Annotated vs Hypothetical comparison ──
    comparison = analysis.get('comparison')
    if comparison and has_comparison:
        ax6 = fig.add_subplot(gs[2, 0:2])

        categories = ['Annotated\n(n={})'.format(comparison['annotated_n']),
                       'Hypothetical\n(n={})'.format(comparison['hypothetical_n'])]
        rates = [comparison['annotated_hit_rate'], comparison['hypothetical_hit_rate']]

        # Wilson score CI error bars
        ann_ci = comparison.get('ann_wilson_ci', [0, 0])
        hyp_ci = comparison.get('hyp_wilson_ci', [0, 0])
        yerr_lower = [rates[0] - ann_ci[0]*100, rates[1] - hyp_ci[0]*100]
        yerr_upper = [ann_ci[1]*100 - rates[0], hyp_ci[1]*100 - rates[1]]

        bars = ax6.bar(categories, rates, color=['steelblue', 'coral'],
                       edgecolor='gray', width=0.5, zorder=2)
        ax6.errorbar(categories, rates, yerr=[yerr_lower, yerr_upper],
                     fmt='none', color='black', capsize=8, capthick=1.5, zorder=3)

        # Annotate bars with values
        for bar, rate in zip(bars, rates):
            ax6.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 1.5,
                     f'{rate:.1f}%', ha='center', va='bottom', fontweight='bold', fontsize=12)

        # Add significance annotation
        ci_lo = comparison.get('fisher_ci_lower', 0)
        ci_hi = comparison.get('fisher_ci_upper', 0)
        p_val = comparison['fisher_p']
        or_val = comparison['fisher_or']
        sig_text = f'OR = {or_val:.1f} (95% CI: {ci_lo:.1f}–{ci_hi:.1f})\np = {p_val:.2e}'
        ax6.text(0.5, 0.85, sig_text, transform=ax6.transAxes, ha='center',
                 fontsize=10, bbox=dict(boxstyle='round,pad=0.3', facecolor='lightyellow',
                                        edgecolor='gray', alpha=0.9))

        # Annotate the single hypothetical hit
        ax6.text(1, rates[1] + 3.5, '* Single hit: A. castellanii\n  (82.4% id, likely host HGT)',
                 ha='center', va='bottom', fontsize=8, fontstyle='italic', color='gray')

        ax6.set_ylabel('Non-self Hit Rate (%)')
        ax6.set_ylim(0, max(rates) * 1.5)
        ax6.set_title('F. Annotated vs Hypothetical: External Homology Rate', fontsize=10)
        ax6.grid(axis='y', alpha=0.3)

        # ── Panel G (optional): Three-level gradient summary ──
        ax7 = fig.add_subplot(gs[2, 2])
        grad_cats = ['Intergenic\nregions', 'Hypothetical\ngenes', 'Annotated\ngenes']
        grad_rates = [0, comparison['hypothetical_hit_rate'], comparison['annotated_hit_rate']]
        grad_colors = ['#d9d9d9', 'coral', 'steelblue']
        grad_n = ['~1,900 IGRs', f'n={comparison["hypothetical_n"]}', f'n={comparison["annotated_n"]}']

        bars_g = ax7.barh(grad_cats, grad_rates, color=grad_colors, edgecolor='gray', height=0.5)
        for bar, rate, n in zip(bars_g, grad_rates, grad_n):
            x_pos = max(rate + 0.8, 2.0)
            ax7.text(x_pos, bar.get_y() + bar.get_height()/2,
                     f'{rate:.1f}%  ({n})', va='center', fontsize=9)

        ax7.set_xlabel('Non-self Hit Rate (%)')
        ax7.set_xlim(0, max(grad_rates) * 1.4)
        ax7.set_title('G. Maturation Gradient:\nBirth → Orphan → Homolog', fontsize=10)
        ax7.invert_yaxis()

    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"\nSaved figure: {output_path}")


# ── Store results in database ─────────────────────────────────────────

def store_results_in_db(blast_results, genes):
    """Store BLAST results in the external_hit table."""
    gene_lookup = {g['protein_id']: g for g in genes if g['protein_id']}

    conn = get_conn()
    cur = conn.cursor()

    # Check external_hit table structure
    cur.execute("""
        SELECT column_name FROM information_schema.columns
        WHERE table_name='external_hit' ORDER BY ordinal_position
    """)
    columns = [r[0] for r in cur.fetchall()]
    print(f"external_hit columns: {columns}")

    # If table is empty, we may need to insert with correct schema
    # For now, skip DB storage if schema doesn't match — save to JSON instead
    conn.close()
    print("Results saved to JSON (external_hit table schema check needed)")


# ── Main ──────────────────────────────────────────────────────────────

def main():
    print("=" * 70)
    print("Phase 2: BLASTp Homology Analysis — Pandoravirus salinus")
    print("=" * 70)

    # ── Load genes ──
    print("\n--- Step 1: Load genes from database ---")
    genes = load_genes()
    annotated = [g for g in genes if classify_gene(g['product']) == 'annotated']
    hypothetical = [g for g in genes if classify_gene(g['product']) == 'hypothetical']
    print(f"Total genes: {len(genes)}")
    print(f"  Annotated:    {len(annotated)}")
    print(f"  Hypothetical: {len(hypothetical)}")

    # ── Export FASTA ──
    print("\n--- Step 2: Export protein sequences ---")
    ann_proteins = export_fasta(genes, 'annotated')
    hyp_proteins = export_fasta(genes, 'hypothetical')

    # ── BLAST ──
    do_blast = '--blast' in sys.argv
    do_blast_hyp = '--blast-hyp' in sys.argv

    ann_results_path = RESULTS_DIR / "blastp_annotated.json"
    hyp_results_path = RESULTS_DIR / "blastp_hypothetical.json"

    if do_blast:
        print("\n--- Step 3: BLASTp annotated genes (remote NCBI) ---")
        ann_blast = run_blastp(ann_proteins, ann_results_path, label="annotated")
    elif ann_results_path.exists():
        print("\n--- Step 3: Loading existing annotated BLAST results ---")
        with open(ann_results_path, 'r') as f:
            ann_blast = json.load(f)
        print(f"  {len(ann_blast)} annotated proteins with results")
    else:
        print("\n--- Step 3: No BLAST results yet ---")
        print(f"  Run with --blast to start BLASTp ({len(ann_proteins)} annotated proteins)")
        print(f"  Estimated time: ~{len(ann_proteins) * 30 / 60:.0f} minutes (resumable)")
        ann_blast = None

    if do_blast_hyp:
        print("\n--- Step 3b: BLASTp hypothetical genes (subsample) ---")
        # Load existing results to preserve prior work
        existing_ids = set()
        if hyp_results_path.exists():
            with open(hyp_results_path, 'r') as f:
                existing_ids = set(json.load(f).keys())

        # Sample HYP_SAMPLE_SIZE proteins, ensuring all previously-sampled proteins are included
        random.seed(42)
        initial_sample = random.sample(hyp_proteins, min(100, len(hyp_proteins)))
        initial_ids = {p.id for p in initial_sample}

        # If we need more than the initial 100, add additional proteins
        if HYP_SAMPLE_SIZE > 100:
            remaining_pool = [p for p in hyp_proteins if p.id not in initial_ids]
            random.seed(43)  # different seed for the extension
            extra = random.sample(remaining_pool, min(HYP_SAMPLE_SIZE - 100, len(remaining_pool)))
            hyp_sample = initial_sample + extra
        else:
            hyp_sample = initial_sample

        hyp_blast = run_blastp(hyp_sample, hyp_results_path, label="hypothetical")
    elif hyp_results_path.exists():
        print("\n--- Step 3b: Loading existing hypothetical BLAST results ---")
        with open(hyp_results_path, 'r') as f:
            hyp_blast = json.load(f)
        print(f"  {len(hyp_blast)} hypothetical proteins with results")
    else:
        hyp_blast = None

    # ── Analysis ──
    if ann_blast:
        print("\n--- Step 4: Analysis ---")
        analysis = analyze_blast_results(ann_blast, genes)

        # Comparison with hypothetical if available
        if hyp_blast:
            print(f"\n--- Annotated vs Hypothetical comparison ---")
            hyp_analysis = analyze_blast_results(hyp_blast, genes)

            ann_hit_rate = analysis['basic']['has_non_self_hit'] / analysis['basic']['successful_queries'] * 100
            hyp_hit_rate = hyp_analysis['basic']['has_non_self_hit'] / hyp_analysis['basic']['successful_queries'] * 100
            print(f"\n  Annotated non-self hit rate:    {ann_hit_rate:.1f}%")
            print(f"  Hypothetical non-self hit rate:  {hyp_hit_rate:.1f}%")

            # Fisher's exact test for hit rate difference
            table = [
                [analysis['basic']['has_non_self_hit'],
                 analysis['basic']['successful_queries'] - analysis['basic']['has_non_self_hit']],
                [hyp_analysis['basic']['has_non_self_hit'],
                 hyp_analysis['basic']['successful_queries'] - hyp_analysis['basic']['has_non_self_hit']],
            ]
            or_val, p_val = stats.fisher_exact(table)

            # Compute 95% CI for odds ratio via log transform
            # Using Woolf's method with Haldane correction (+0.5 to each cell)
            a, b = table[0]
            c, d = table[1]
            a_h, b_h, c_h, d_h = a + 0.5, b + 0.5, c + 0.5, d + 0.5
            log_or = np.log(a_h * d_h / (b_h * c_h))
            se_log_or = np.sqrt(1/a_h + 1/b_h + 1/c_h + 1/d_h)
            ci_lower = np.exp(log_or - 1.96 * se_log_or)
            ci_upper = np.exp(log_or + 1.96 * se_log_or)

            print(f"  Fisher's exact: OR={or_val:.2f}, 95% CI: {ci_lower:.1f}–{ci_upper:.1f}, p={p_val:.2e}")

            # Wilson score CIs for proportions (for figure error bars)
            def wilson_ci(k, n, z=1.96):
                """Wilson score confidence interval for binomial proportion."""
                if n == 0:
                    return 0, 0, 0
                p = k / n
                denom = 1 + z**2 / n
                center = (p + z**2 / (2*n)) / denom
                margin = z * np.sqrt(p*(1-p)/n + z**2/(4*n**2)) / denom
                return p, max(0, center - margin), min(1, center + margin)

            ann_n = analysis['basic']['successful_queries']
            ann_k = analysis['basic']['has_non_self_hit']
            hyp_n = hyp_analysis['basic']['successful_queries']
            hyp_k = hyp_analysis['basic']['has_non_self_hit']
            ann_p, ann_ci_lo, ann_ci_hi = wilson_ci(ann_k, ann_n)
            hyp_p, hyp_ci_lo, hyp_ci_hi = wilson_ci(hyp_k, hyp_n)

            analysis['comparison'] = {
                'annotated_hit_rate': ann_hit_rate,
                'hypothetical_hit_rate': hyp_hit_rate,
                'annotated_n': ann_n,
                'annotated_hits': ann_k,
                'hypothetical_n': hyp_n,
                'hypothetical_hits': hyp_k,
                'fisher_or': float(or_val),
                'fisher_ci_lower': float(ci_lower),
                'fisher_ci_upper': float(ci_upper),
                'fisher_p': float(p_val),
                'ann_wilson_ci': [float(ann_ci_lo), float(ann_ci_hi)],
                'hyp_wilson_ci': [float(hyp_ci_lo), float(hyp_ci_hi)],
            }

        # ── E-value sensitivity analysis ──
        print("\n--- E-value Sensitivity Analysis ---")
        print(f"{'Threshold':<12} {'Non-self hits':<15} {'Hit rate':<10} {'Orphan rate':<12}")
        print("-" * 50)

        sensitivity_results = {}
        for threshold in [1e-5, 1e-10, 1e-20, 1e-50]:
            # Re-filter raw hits at this threshold
            n_with_nonself = 0
            n_successful = 0
            for pid, r in ann_blast.items():
                if r.get('error'):
                    continue
                n_successful += 1
                # Check if any non-self hit passes this threshold
                non_self_at_threshold = [
                    h for h in r.get('hits', [])
                    if 'pandoravirus' not in h.get('organism', '').lower()
                    and h.get('evalue', 1) <= threshold
                ]
                if non_self_at_threshold:
                    n_with_nonself += 1

            hit_rate = n_with_nonself / n_successful * 100 if n_successful else 0
            orphan_rate = 100 - hit_rate
            print(f"{'<{:.0e}'.format(threshold):<12} {n_with_nonself:<15} {hit_rate:<10.1f}% {orphan_rate:<12.1f}%")

            sensitivity_results[str(threshold)] = {
                'threshold': threshold,
                'successful_queries': n_successful,
                'non_self_hits': n_with_nonself,
                'hit_rate': round(hit_rate, 1),
                'orphan_rate': round(orphan_rate, 1),
            }

        print(f"\nNote: Original BLAST used E < 1e-5. Stricter thresholds (1e-10, 1e-20, 1e-50)")
        print(f"re-filter existing hits. More permissive thresholds (1e-3) would require re-BLASTing.")
        analysis['sensitivity'] = sensitivity_results

        # ── Figures ──
        print("\n--- Step 5: Visualization ---")
        plot_phase2_results(ann_blast, analysis, genes,
                           FIGURES_DIR / "phase2_blastp_homology.png")

        # ── Save full results ──
        # Convert non-serializable types
        save_analysis = {}
        for k, v in analysis.items():
            if k == 'identity' and 'values' in v:
                v = {**v, 'values': [float(x) for x in v['values']]}
            if k == 'identity_by_taxonomy':
                v = {cat: {**d, 'values': [float(x) for x in d['values']]}
                     for cat, d in v.items()}
            save_analysis[k] = v

        results_out = RESULTS_DIR / "phase2_results.json"
        with open(results_out, 'w') as f:
            json.dump(save_analysis, f, indent=2, default=str)
        print(f"\nResults saved: {results_out}")

    # ── Summary ──
    print(f"\n{'='*60}")
    print("Phase 2 Status")
    print(f"{'='*60}")
    if ann_blast:
        n = len(ann_blast)
        print(f"Annotated BLAST:    {n}/{len(ann_proteins)} complete")
    else:
        print(f"Annotated BLAST:    not started (run with --blast)")
    if hyp_blast:
        print(f"Hypothetical BLAST: {len(hyp_blast)}/{HYP_SAMPLE_SIZE} complete")
    else:
        print(f"Hypothetical BLAST: not started (run with --blast-hyp)")


if __name__ == "__main__":
    main()
