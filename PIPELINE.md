# Pandoravirus 2.0 — Reproducible Analysis Pipeline
## Audit-Ready Module Specification

---

## Overview

This document specifies the complete analytical pipeline for "Computational evidence
for a distributed AT-based regulatory mechanism enabling de novo gene birth in
Pandoravirus salinus." Each module is defined with explicit inputs, outputs, parameters,
verification checks, and dependencies, enabling reproduction by human researchers or
automated agents.

**Pipeline DAG (Directed Acyclic Graph):**

```
[GenBank Accession] ──► MODULE 1 (Data Acquisition)
        │
        ▼
[PostgreSQL Database] ──► MODULE 2 (Genome Architecture)
        │                         │
        ▼                         ▼
MODULE 3 (Null Models) ──► MODULE 4 (Codon Usage)
        │                         │
        ▼                         ▼
MODULE 5 (Gene Prediction) ──► MODULE 6 (Evolutionary Classification)
        │                              │
        ▼                              ▼
MODULE 7 (Gene Families + Proto-genes) ──► MODULE 8 (Protein Structure)
        │
        ▼
MODULE 9 (Regulatory Signals) ──► MODULE 10 (Transcriptomic Validation)
        │                                    │
        ▼                                    ▼
MODULE 11 (BLASTp Homology) ──► MODULE 12 (Integration + Figures)
```

**Global Requirements:**
- Python 3.11+
- PostgreSQL 15+
- BioPython 1.86+
- NumPy, SciPy, scikit-learn, matplotlib, seaborn, pandas
- Internet access (for NCBI BLAST and Entrez queries in Modules 1 and 11)

---

## MODULE 1: Data Acquisition & Database Setup

**Purpose:** Obtain genome sequences and annotations, build the analytical database.

### Inputs
| Input | Source | Identifier |
|-------|--------|------------|
| P. salinus genome + annotations | NCBI GenBank | NC_022098.1 |
| PhiX174 genome | NCBI GenBank | NC_001422.1 |
| Lambda genome | NCBI GenBank | NC_001416.1 |
| Mimivirus genome | NCBI GenBank | NC_014649.1 |
| P. salinus GFF (RNA-seq informed) | IGS GBrowse | http://www.igs.cnrs-mrs.fr/pandoraviruses/ |
| A. castellanii codon table | Kazusa | Species ID 5755 |

### Process
1. Download genome sequences and annotations via Entrez
2. Parse GenBank flat files: extract CDS, gene coordinates, translations, products
3. Create PostgreSQL database with schema (22 tables)
4. Load genome sequences as position-indexed nucleotide arrays
5. Load gene annotations with coordinates, strand, product, protein_id
6. Compute per-gene: GC content, gene length, intergenic region coordinates
7. Download and parse GFF from IGS GBrowse for Phase 10

### Outputs
| Output | Format | Description | Verification |
|--------|--------|-------------|-------------|
| `pandoravirus.sql` | SQL dump | Complete database | 22 tables, 1,430 CDS for P. salinus |
| `data/genomes/` | FASTA + GBK | Raw genome files | 4 genomes, checksums match NCBI |
| `data/gff/` | GFF3 | RNA-seq-informed annotations | Gene count matches Legendre 2018 |
| `data/codon_tables/` | JSON | A. castellanii codon frequencies | 64 codons, frequencies sum to 1.0 |

### Verification Checks
- [ ] P. salinus genome length = 2,473,870 bp
- [ ] P. salinus CDS count = 1,430
- [ ] GC content = 61.7% (±0.1%)
- [ ] All 4 comparison genomes loaded
- [ ] Database integrity: foreign key constraints pass

### Parameters
| Parameter | Value | Justification |
|-----------|-------|---------------|
| Database engine | PostgreSQL 15 | Spatial indexing for position queries |
| Genome ID (P. salinus) | 4 | Internal database key |

---

## MODULE 2: Genome Architecture Analysis

**Purpose:** Characterize basic compositional features across all four genomes.

### Inputs
| Input | Source |
|-------|--------|
| Populated PostgreSQL database | Module 1 |

### Process
1. Compute 8 genome-wide metrics: GC content, coding density, strand bias, intergenic length distribution, Shannon entropy, gene length distribution, sliding-window GC (1kb windows, 100bp step)
2. Compute per-gene metrics: GC, GC3, GC12, gene length, upstream/downstream intergenic lengths
3. Compare all metrics across 4-genome ladder
4. Calculate coefficient of variation for GC across sliding windows

### Outputs
| Output | Format | Description | Verification |
|--------|--------|-------------|-------------|
| `results/phase2_architecture.json` | JSON | All computed metrics | P. salinus GC CV = 6.1% |
| `figures/genome_architecture.png` | PNG | Comparative panel figure | 4 genomes shown |

### Verification Checks
- [ ] P. salinus coding density ~80%
- [ ] P. salinus GC CV = 6.1% (lowest of 4 genomes)
- [ ] Intergenic region count ~1,900
- [ ] Sliding window GC range reasonable (no extreme outliers from assembly errors)

---

## MODULE 3: Null Model Calibration

**Purpose:** Build Markov-1 null models for k-mer analysis; establish baseline expectations.

### Inputs
| Input | Source |
|-------|--------|
| Genome sequences | Module 1 (database) |

### Process
1. Compute observed k-mer frequencies (k = 4, 6, 8) for each genome
2. Generate Markov-1 expected frequencies preserving dinucleotide composition
3. Compute z-scores: z = (observed - expected) / std_expected
4. Apply dual filtering: z-score significance (|z| > 3) AND effect-size (O/E < 0.5 or > 2.0)
5. Report: fraction of z > 3 signals that lack meaningful effect sizes

### Outputs
| Output | Format | Description | Verification |
|--------|--------|-------------|-------------|
| `results/phase3_null_models.json` | JSON | k-mer O/E ratios, z-scores | AGCT 1,500-fold depletion confirmed |
| `data/null_models/markov1_*.json` | JSON | Per-genome null frequency tables | Frequencies sum to 1.0 per k |

### Verification Checks
- [ ] AGCT tetranucleotide O/E ratio < 0.001 in P. salinus (confirms Poirot et al. 2019)
- [ ] 88.6% of z > 3 signals lack effect size > 2.0 (validates dual-filter approach)
- [ ] Null model dinucleotide frequencies match observed (by construction)

### Parameters
| Parameter | Value | Justification |
|-----------|-------|---------------|
| k values | 4, 6, 8 | Standard range for regulatory motif detection |
| z-score threshold | 3.0 | Conventional significance |
| Effect-size threshold | 2.0 (O/E) | Eliminates statistically significant but biologically trivial signals |

---

## MODULE 4: Codon Usage Analysis

**Purpose:** Characterize per-gene codon usage patterns; compare annotated vs. ORFan genes.

### Inputs
| Input | Source |
|-------|--------|
| Gene translations and CDS sequences | Module 1 (database) |
| A. castellanii codon frequencies | Module 1 (Kazusa) |

### Process
1. Compute per-gene: RSCU (64 values), GC3, GC12, ENC, CAI
2. CAI reference: A. castellanii codon frequencies from Kazusa
3. Compare annotated vs. ORFan genes (Wilcoxon rank-sum tests)
4. Compute codon usage distance metrics for proto-gene disambiguation (Module 7)

### Outputs
| Output | Format | Description | Verification |
|--------|--------|-------------|-------------|
| `results/phase4_codon_usage.json` | JSON | Per-gene RSCU, GC3, GC12, ENC, CAI | 1,430 genes with complete metrics |
| Database table `codon_usage` | SQL | Indexed per-gene codon statistics | Queryable by gene_id |

### Verification Checks
- [ ] ORFan GC3 difference from annotated = +0.015 (not significant → compositional homogeneity)
- [ ] ORFan CAI difference from annotated = +0.008 (not significant)
- [ ] All CAI values in range [0, 1]
- [ ] ENC values in range [20, 61]

---

## MODULE 5: Gene Prediction Validation

**Purpose:** Compare GenBank annotations against ab initio predictions; identify novel candidate ORFs.

### Inputs
| Input | Source |
|-------|--------|
| P. salinus genome sequence | Module 1 |
| GenBank CDS annotations | Module 1 |

### Process
1. Run Prodigal (v2.6.3, meta mode) on P. salinus genome
2. Compare Prodigal predictions vs. GenBank annotations
3. Match criterion: ≥80% reciprocal overlap
4. Classify: confirmed (both), GenBank-only, Prodigal-only
5. Characterize Prodigal-only predictions: length, GC, position

### Outputs
| Output | Format | Description | Verification |
|--------|--------|-------------|-------------|
| `data/prodigal/` | GFF + FAA | Raw Prodigal output | Reproducible with same Prodigal version |
| `results/phase5_prediction.json` | JSON | Classification counts, Prodigal-only characteristics | 347 novel Prodigal-only predictions |

### Verification Checks
- [ ] Prodigal total predictions > GenBank count (expected: Prodigal predicts more)
- [ ] 347 Prodigal-only predictions identified
- [ ] Prodigal-only mean GC < GenBank mean GC (expected: less mature)
- [ ] Prodigal-only mean length < GenBank mean length

### Parameters
| Parameter | Value | Justification |
|-----------|-------|---------------|
| Prodigal mode | meta | Single genome, no training set |
| Overlap threshold | 80% reciprocal | Standard for gene prediction comparison |

---

## MODULE 6: Evolutionary Origin Classification

**Purpose:** Cluster genes by compositional features; define the native zone.

### Inputs
| Input | Source |
|-------|--------|
| Per-gene GC, GC3, GC12, ENC, CAI | Module 4 |
| Gene positions | Module 1 |

### Process
1. Feature matrix: 5 features × n genes (GC, GC3, GC12, ENC, CAI)
2. Gaussian Mixture Model clustering (BIC-selected component count)
3. Assess spatial distribution of clusters
4. Define native zone boundaries from cluster enrichment pattern
5. Compare cluster composition to identify native vs. foreign-like genes

### Outputs
| Output | Format | Description | Verification |
|--------|--------|-------------|-------------|
| `results/phase6_classification.json` | JSON | Cluster assignments, native zone boundaries | Native zone = 1.3-2.0 Mb |

### Verification Checks
- [ ] GMM identifies compositionally distinct clusters
- [ ] Native zone (1.3-2.0 Mb) shows enrichment for a specific cluster
- [ ] Zone boundaries consistent with Christo-Foroux/Aherfi pangenome observations

### Parameters
| Parameter | Value | Justification |
|-----------|-------|---------------|
| GMM components | BIC-selected | Data-driven model selection |
| Native zone | 1.3-2.0 Mb | Derived from clustering + literature |

---

## MODULE 7: Gene Families & Proto-Gene Classification

**Purpose:** Identify gene families; classify all ORFs into the five-tier proto-gene continuum.

### Inputs
| Input | Source |
|-------|--------|
| All protein sequences (GenBank + Prodigal) | Modules 1, 5 |
| Codon usage metrics | Module 4 |
| Gene prediction classifications | Module 5 |

### Process
1. All-vs-all protein 3-mer frequency vectors (1,021,735 pairwise comparisons)
2. Cosine similarity matrix
3. Hierarchical clustering to identify gene families
4. Null calibration: composition-matched random proteins (confirm diffuse cloud is artifact)
5. Classify into 5 tiers based on prediction agreement:
   - T1: Annotated, Prodigal-confirmed (n=516)
   - T2: Hypothetical, Prodigal-confirmed (n=784)
   - T3: GenBank-only, Prodigal missed (n=130)
   - T4: Prodigal high-confidence only (n=335)
   - T5: Prodigal low-confidence only (n=527)
6. Proto-gene codon disambiguation: compare T4/T5 codon usage vs. length-matched random intergenic ORFs
7. Proto-gene spatial analysis: test for clustering near established genes

### Outputs
| Output | Format | Description | Verification |
|--------|--------|-------------|-------------|
| `results/phase7_families.json` | JSON | 14 gene families, similarity statistics | Mean similarity 0.109, null 0.105 |
| `results/phase7_protogenes.json` | JSON | 5-tier classification, codon tests | T5 codon p = 0.0005 |
| `results/phase7_continuum.json` | JSON | Tier-level summary statistics (Table 1) | Monotonic GC/length gradient |

### Verification Checks
- [ ] 14 gene families detected (diffuse cloud confirmed as artifact: Cohen's d = 0.06)
- [ ] T1 through T5 show monotonic decline in GC (0.650 → 0.547)
- [ ] T1 through T5 show monotonic decline in length (1,452 → 216 bp)
- [ ] T4 codon usage vs. random: p = 0.020
- [ ] T5 codon usage vs. random: p = 0.0005
- [ ] T4 spatial clustering: p = 0.99 (no clustering)
- [ ] T5 spatial clustering: p = 1.0 (no clustering)
- [ ] Spearman ρ = -1.0 for GC across tiers

---

## MODULE 8: Protein Structure Prediction

**Purpose:** Assess structural relationships between proto-genes and established proteins.

### Inputs
| Input | Source |
|-------|--------|
| Representative protein sequences (170 proteins spanning all tiers) | Module 7 |

### Process
1. Generate ESM-2 (650M parameter) embeddings for 170 representative proteins
2. PCA reduction: top 20 components (target: ~85% variance explained)
3. Hierarchical clustering on PCA-reduced embeddings
4. Null calibration: confirm structural signal is genuine
5. Assess GC-structure correlation (PC1 vs. GC content)
6. Identify structural nearest-neighbors for proto-gene candidates

### Outputs
| Output | Format | Description | Verification |
|--------|--------|-------------|-------------|
| `results/phase8_structure.json` | JSON | Clustering results, null comparison | Within-family sim 0.533 vs null 0.368 |
| `data/embeddings/esm2_*.npy` | NumPy | Raw ESM-2 embeddings | 170 proteins × embedding dimension |

### Verification Checks
- [ ] PCA top 20 components explain ~85.8% variance
- [ ] Within-family structural similarity (0.533) >> null (0.368); Cohen's d = 1.45
- [ ] 41% of PC1 variance correlates with GC content (known confound, reported)
- [ ] Proto-gene nearest-neighbor structural similarity p = 0.90 vs established genes

### Parameters
| Parameter | Value | Justification |
|-----------|-------|---------------|
| ESM-2 model | 650M parameters | Balance of quality and compute |
| PCA components | 20 | Capture ~85% variance |
| Representative sample | 170 proteins | Spans all tiers and clusters |

---

## MODULE 9: Regulatory Signal Discovery

**Purpose:** Characterize the AT-based regulatory landscape of P. salinus.

### Inputs
| Input | Source |
|-------|--------|
| Genome sequence (position-indexed) | Module 1 |
| Gene annotations (coordinates, strand) | Module 1 |
| Intergenic region coordinates | Module 2 |
| Proto-gene classifications | Module 7 |

### Process
1. Inventory all poly-A/T runs ≥ 5 bp genome-wide
2. Classify runs as coding/intergenic
3. Meta-gene AT content profiles: ±500 bp around all gene starts/ends
4. Upstream k-mer enrichment analysis
5. Intergenic grammar: compare convergent vs. divergent regions (length, run density)
6. Strand-asymmetric boundary positioning analysis
7. ORFan vs. annotated regulatory signature comparison
8. Proto-gene regulatory context: compare T4/T5 upstream AT to annotated gene upstream AT

### Outputs
| Output | Format | Description | Verification |
|--------|--------|-------------|-------------|
| `results/phase9_regulatory.json` | JSON | All regulatory signal metrics | Start ΔAT = +7.0 pp, End ΔAT = +9.5 pp |
| `results/phase9_polyAT_runs.csv` | CSV | Position, length, strand, type for every run | ~14,000 runs total |
| `figures/regulatory_signals.png` | PNG | Meta-gene profiles, grammar comparison | |

### Verification Checks
- [ ] Poly-A runs: ~6,955; Poly-T runs: ~7,010
- [ ] Intergenic fraction: 82.8–83.9%
- [ ] Meta-gene AT spike at starts: +7.0 pp (unfiltered)
- [ ] Meta-gene AT spike at ends: +9.5 pp (unfiltered)
- [ ] Divergent regions 17% longer than convergent (p < 0.0001)
- [ ] T4 upstream AT: 45.0% (p = 0.14 vs annotated) — not significantly different
- [ ] T5 upstream AT: 46.7% (p = 0.25 vs annotated) — not significantly different
- [ ] P. salinus boundary signal >> Lambda (+1 pp) and Mimivirus (+1-2 pp)

---

## MODULE 10: Transcriptomic Validation

**Purpose:** Validate computationally predicted regulatory signals against experimentally informed annotations.

### Inputs
| Input | Source |
|-------|--------|
| Intergenic poly-A/T run inventory | Module 9 |
| RNA-seq-informed GFF annotations | Module 1 (IGS GBrowse) |
| Gene coordinates with strand information | Module 1 |

### Process
1. Filter poly-A/T runs to intergenic only (remove 2,325 coding-overlapping runs)
2. Strand-aware boundary classification:
   - Promoter: poly-A upstream of sense-strand start, OR poly-T upstream of antisense-strand start
   - Terminator: poly-T downstream of sense-strand end, OR poly-A downstream of antisense-strand end
3. Compute fraction in correct regulatory position
4. Binomial test against 50% null
5. Recompute meta-gene AT profiles on filtered intergenic runs
6. Hairpin null model: strict parameters (stem ≥ 6 bp, exact complementarity)
   - Scan 1,430 downstream regions for palindromic sequences
   - Compare observed rate to Markov-1 null expectation (1,000 permutations)
7. Two-population test: compare boundary-proximal (≤200 bp) vs. dispersed runs
   - Wilcoxon rank-sum on AT ratio and length

### Outputs
| Output | Format | Description | Verification |
|--------|--------|-------------|-------------|
| `results/phase1b_results.json` | JSON | All validation metrics | 74% correct positioning |
| `figures/phase1b_corrected.png` | PNG | 6-panel validation figure | |

### Verification Checks
- [ ] Intergenic runs after filtering: 11,640 (from ~13,965)
- [ ] Coding runs removed: 2,325 (16.6%)
- [ ] Strand-aware correct positioning: 74% (binomial p < 10⁻⁶)
- [ ] Promoter-consistent runs: 2,332
- [ ] Terminator-consistent runs: 1,227
- [ ] Promoter:terminator ratio: 1.9:1
- [ ] Meta-gene AT spike (filtered): +10.3 pp at starts, +9.8 pp at ends
- [ ] Hairpin observed rate: 11.6%
- [ ] Hairpin null rate: 9.7%
- [ ] Hairpin excess: 2.0 pp (z = 3.6)
- [ ] Two-population AT ratio test: p = 0.105 (not significant → identical populations)
- [ ] Two-population length test: p = 0.093 (not significant → identical populations)

### Parameters
| Parameter | Value | Justification |
|-----------|-------|---------------|
| Intergenic filter | Remove runs overlapping any CDS | Conservative — avoids counting coding poly-A/T |
| Boundary proximity | 200 bp | Captures regulatory-range signals |
| Hairpin stem minimum | 6 bp | Strict — avoids short palindrome noise |
| Hairpin complementarity | Exact | Conservative — no mismatches allowed |
| Null permutations | 1,000 | Sufficient for z-score stability |

---

## MODULE 11: BLASTp Homology Analysis

**Purpose:** Characterize external homology landscape; test maturation model predictions.

### Inputs
| Input | Source |
|-------|--------|
| Annotated protein sequences (528) | Module 1 (database) |
| Hypothetical protein sequences (902 total) | Module 1 (database) |
| Gene positions and product annotations | Module 1 |
| Native zone boundaries | Module 6 |

### Process
1. Export annotated protein FASTA (528 proteins, non-hypothetical)
2. Remote BLASTp against NCBI nr database (resumable, one query at a time)
3. Parse XML results: extract top 10 hits per query
4. Filter non-self hits (exclude Pandoravirus strains)
5. Taxonomic classification:
   a. Curated organism-level lookup table (50 unique organisms)
   b. Ordered keyword fallback (giant viruses checked before generic "virus")
   c. NCBI Entrez API lookup for unrecognized organisms
6. Export hypothetical protein FASTA subsample:
   - Initial 100 (random seed 42) + additional 200 (random seed 43) = 300 total
7. BLASTp hypothetical subsample under identical parameters
8. Compute: identity distribution, spatial enrichment (Fisher's exact), product-level hit rates
9. Annotated vs. hypothetical comparison (Fisher's exact + Haldane-corrected OR CI)
10. E-value sensitivity: re-filter existing results at 10⁻¹⁰, 10⁻²⁰, 10⁻⁵⁰

### Outputs
| Output | Format | Description | Verification |
|--------|--------|-------------|-------------|
| `results/blastp_annotated.json` | JSON | Per-protein BLAST results (528 entries) | 521 successful, 7 errors |
| `results/blastp_hypothetical.json` | JSON | Per-protein BLAST results (300 entries) | 299 successful, 1 error |
| `results/phase2_results.json` | JSON | All analysis metrics | See verification checks |
| `figures/phase2_blastp_homology.png` | PNG | 5+ panel figure | |

### Verification Checks

**Basic statistics:**
- [ ] Annotated queries attempted: 528
- [ ] Annotated successful: 521 (7 connection errors excluded)
- [ ] Annotated non-self hits: 95 (18.2%)
- [ ] Annotated self-only: 426 (81.8%)
- [ ] Genuine no-hits: 0
- [ ] Hypothetical queries attempted: 300
- [ ] Hypothetical successful: 299 (1 error)
- [ ] Hypothetical non-self hits: 9 (3.0%)

**Taxonomy (n=95 annotated non-self hits):**
- [ ] Marine algae: 28 (29.5%)
- [ ] Giant virus: 24 (25.3%)
- [ ] Amoebozoa: 13 (13.7%)
- [ ] Bacteria: 11 (11.6%)
- [ ] Insect: 11 (11.6%)
- [ ] Other eukaryote: 7 (7.4%)
- [ ] Other virus: 1 (1.1%)
- [ ] No prasinophytes classified as bacteria (taxonomy correction verified)

**Identity distribution:**
- [ ] Mean: 44.7%
- [ ] Median: 42.3%
- [ ] Range: 25.8%–84.0%

**Spatial enrichment:**
- [ ] Native zone hits: 45/95
- [ ] Fisher's exact OR: 1.94
- [ ] Fisher's exact p: 0.004

**Annotated vs. hypothetical comparison:**
- [ ] Fisher's exact OR: 7.19
- [ ] 95% CI: 3.5–13.6
- [ ] p-value: 7.79 × 10⁻¹²

**E-value sensitivity:**
- [ ] 1e-5: 95 hits (18.2% non-self rate, 81.8% orphan rate)
- [ ] 1e-10: 88 hits (16.9%, 83.1%)
- [ ] 1e-20: 83 hits (15.9%, 84.1%)
- [ ] 1e-50: 50 hits (9.6%, 90.4%)

**Hypothetical hit taxonomy (n=9):**
- [ ] Giant virus: 6 (66.7%)
- [ ] Amoebozoa: 2 (22.2%)
- [ ] Bacteria: 1 (11.1%)
- [ ] Marine algae: 0
- [ ] Insect: 0
- [ ] 7 of 9 in native zone (OR=11.9, p=0.0009)

### Parameters
| Parameter | Value | Justification |
|-----------|-------|---------------|
| BLAST database | nr (NCBI non-redundant protein) | Most comprehensive protein database |
| E-value threshold | 1e-5 (primary) | Standard for distant homology; sensitivity tested at 1e-10, 1e-20, 1e-50 |
| Max hits per query | 10 | Sufficient for taxonomy; top hit used for primary analysis |
| BLAST method | Remote NCBI qblast | Ensures current database; no local installation needed |
| Query delay | 1 second between queries | Respects NCBI rate limits |
| Random seed (hyp initial) | 42 | Reproducible subsample |
| Random seed (hyp expansion) | 43 | Reproducible expansion |
| Non-self filter | Exclude "pandoravirus" in organism | Removes all Pandoravirus strains |

### Known Issues & Mitigations
| Issue | Mitigation |
|-------|-----------|
| NCBI nr database updates over time | Report database access date in manuscript; results may vary slightly with newer database versions |
| Connection failures (7 annotated, 1 hypothetical) | Excluded from denominators; reported as connection errors, not genuine no-hits |
| Taxonomy misclassification risk | Curated lookup table with 50 organisms; all classifications manually verified |
| Remote BLAST non-determinism | Results are resumable/cached; re-running may produce slightly different E-values due to database growth |

---

## MODULE 12: Integration & Figure Generation

**Purpose:** Combine all module outputs into final manuscript figures and tables.

### Inputs
| Input | Source |
|-------|--------|
| Proto-gene continuum data | Module 7 |
| Regulatory signal metrics | Module 9 |
| Transcriptomic validation results | Module 10 |
| BLASTp homology results | Module 11 |
| All intermediate results | Modules 2-8 |

### Process
1. Generate Figure 1: Proto-gene continuum (4 panels)
2. Generate Figure 2: AT regulatory system + validation (6 panels)
3. Generate Figure 3: BLASTp homology landscape (5 panels)
4. Generate Figure 4: Maturation gradient (3 panels)
5. Generate Tables 1-4 for manuscript
6. Generate supplementary figures S1, S2
7. Generate supplementary tables S1-S3
8. Compile all headline statistics for abstract/results cross-check

### Outputs
| Output | Format | Description |
|--------|--------|-------------|
| `figures/figure1_continuum.png` | PNG (300 dpi) | 4-panel proto-gene continuum |
| `figures/figure2_regulatory.png` | PNG (300 dpi) | 6-panel AT regulatory system |
| `figures/figure3_homology.png` | PNG (300 dpi) | 5-panel BLASTp landscape |
| `figures/figure4_maturation.png` | PNG (300 dpi) | 3-panel maturation gradient |
| `figures/figS1_phase2_full.png` | PNG (150 dpi) | 7-panel supplementary |
| `figures/figS2_phase1b_full.png` | PNG (150 dpi) | 6-panel supplementary |
| `manuscript/tables.md` | Markdown | Tables 1-4 + S1-S3 |
| `manuscript/statistics_audit.json` | JSON | Every statistic cited in manuscript with source module |

### Verification: Statistics Audit Trail
The `statistics_audit.json` file maps every number cited in the manuscript to its source:

```json
{
  "abstract": {
    "orfan_fraction": {
      "value": "92.5%",
      "source": "module1",
      "computation": "1323_orfans / 1430_total_genes",
      "file": "results/phase2_architecture.json"
    },
    "annotated_orphan_rate": {
      "value": "81.8%",
      "source": "module11",
      "computation": "426 / 521",
      "file": "results/phase2_results.json",
      "path": "basic.self_only / basic.successful_queries"
    },
    "hypothetical_orphan_rate": {
      "value": "97.0%",
      "source": "module11",
      "computation": "290 / 299",
      "file": "results/phase2_results.json",
      "path": "comparison.hypothetical_hit_rate → 100 - 3.0"
    },
    "fisher_or": {
      "value": "7.19",
      "source": "module11",
      "computation": "fisher_exact([[95,426],[9,290]])",
      "file": "results/phase2_results.json",
      "path": "comparison.fisher_or"
    }
  }
}
```

This file enables automated verification: an agent can parse the manuscript, extract every cited statistic, look up its source in the audit trail, recompute from the source data, and flag any discrepancies.

---

## REPOSITORY STRUCTURE

```
pandoravirus-2.0/
├── README.md                          # Project overview + quickstart
├── PIPELINE.md                        # This document
├── LICENSE                            # MIT or CC-BY-4.0
├── environment.yml                    # Conda environment specification
├── requirements.txt                   # pip requirements (fallback)
│
├── data/
│   ├── genomes/                       # Raw genome files (Module 1)
│   ├── gff/                           # RNA-seq-informed annotations (Module 1)
│   ├── codon_tables/                  # Reference codon frequencies (Module 1)
│   ├── prodigal/                      # Ab initio predictions (Module 5)
│   ├── embeddings/                    # ESM-2 embeddings (Module 8)
│   ├── null_models/                   # Markov-1 null tables (Module 3)
│   └── proteins/                      # Exported FASTA files (Module 11)
│
├── scripts/
│   ├── module01_data_acquisition.py
│   ├── module02_genome_architecture.py
│   ├── module03_null_models.py
│   ├── module04_codon_usage.py
│   ├── module05_gene_prediction.py
│   ├── module06_evolutionary_classification.py
│   ├── module07_gene_families.py
│   ├── module08_protein_structure.py
│   ├── module09_regulatory_signals.py
│   ├── module10_transcriptomic_validation.py
│   ├── module11_blastp_homology.py
│   ├── module12_integration_figures.py
│   └── verify_pipeline.py            # Runs all verification checks
│
├── results/
│   ├── module_outputs/                # One JSON per module
│   └── statistics_audit.json          # Manuscript ↔ data mapping
│
├── figures/
│   ├── figure1_continuum.png
│   ├── figure2_regulatory.png
│   ├── figure3_homology.png
│   ├── figure4_maturation.png
│   └── supplementary/
│
├── manuscript/
│   ├── paper_draft.md
│   ├── tables.md
│   └── references.bib
│
└── tests/
    ├── test_module01.py               # Unit tests per module
    ├── ...
    └── test_integration.py            # End-to-end verification
```

---

## EXECUTION GUIDE

### For Human Researchers
```bash
# 1. Set up environment
conda env create -f environment.yml
conda activate pandoravirus

# 2. Initialize database
createdb pandoravirus
python scripts/module01_data_acquisition.py

# 3. Run pipeline (sequential, ~4-6 hours excluding BLAST)
for i in $(seq 2 12); do
    python scripts/module$(printf "%02d" $i)_*.py
done

# 4. Run BLAST (requires internet, ~24-48 hours)
python scripts/module11_blastp_homology.py --blast
python scripts/module11_blastp_homology.py --blast-hyp

# 5. Verify all results
python scripts/verify_pipeline.py
```

### For AI Agents
```yaml
pipeline:
  name: pandoravirus-2.0
  description: De novo gene birth analysis of P. salinus
  
  prerequisites:
    - python: ">=3.11"
    - postgresql: ">=15"
    - internet: true  # Required for BLAST (Module 11)
    
  execution_order:
    - module: 1
      script: module01_data_acquisition.py
      requires_internet: true
      expected_runtime: "5 minutes"
      verify: "SELECT COUNT(*) FROM gene WHERE genome_id=4 -- expect 1430"
      
    - module: 2
      script: module02_genome_architecture.py
      depends_on: [1]
      expected_runtime: "2 minutes"
      verify: "results/phase2_architecture.json exists AND gc_cv == 0.061"
      
    - module: 3
      script: module03_null_models.py
      depends_on: [1]
      expected_runtime: "5 minutes"
      verify: "AGCT O/E ratio < 0.001"
      
    # ... (modules 4-10 follow same pattern)
    
    - module: 11
      script: module11_blastp_homology.py --blast --blast-hyp
      depends_on: [1, 6]
      expected_runtime: "24-48 hours"
      requires_internet: true
      resumable: true
      verify: |
        basic.successful_queries == 521
        basic.has_non_self_hit == 95
        comparison.fisher_or BETWEEN 6.0 AND 9.0
        comparison.fisher_p < 1e-10
      note: "BLAST results may vary slightly with database updates. OR and p-value should remain in similar range."
      
    - module: 12
      script: module12_integration_figures.py
      depends_on: [2,3,4,5,6,7,8,9,10,11]
      expected_runtime: "2 minutes"
      verify: "All 4 main figures + statistics_audit.json generated"
      
  final_verification:
    script: verify_pipeline.py
    description: "Runs all module verification checks and cross-references statistics_audit.json against manuscript"
```

---

## VERSIONING & AUDIT LOG

Each pipeline execution should generate an audit log:

```json
{
  "execution_date": "2026-03-15",
  "pipeline_version": "2.0.0",
  "environment": {
    "python": "3.11.x",
    "biopython": "1.86",
    "scipy": "1.12.x",
    "postgresql": "15.x",
    "prodigal": "2.6.3",
    "esm2_model": "esm2_t33_650M_UR50D"
  },
  "ncbi_blast_database_date": "2026-03-01",
  "modules_executed": [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12],
  "verification_checks_passed": 87,
  "verification_checks_failed": 0,
  "total_runtime_hours": 28.5,
  "notes": "Initial full pipeline execution"
}
```

---

## NOTES ON REPRODUCIBILITY EXPECTATIONS

BLAST results are inherently non-reproducible across time because the nr database grows
continuously. A BLAST search run in 2026 may find hits that did not exist in 2025, or may
return different E-values for the same alignments due to database size changes affecting
the statistical model. The key results (orphan rates, OR values, p-values) should be
robust to small database changes, but exact hit counts may differ.

All other modules (1-10, 12) are fully deterministic given the same input genome sequence
and the same software versions. The random seeds for the hypothetical subsample (42, 43)
ensure that the same proteins are selected for BLAST comparison.

The verification checks include tolerance ranges where appropriate (e.g., Module 11
verifies that the OR is "BETWEEN 6.0 AND 9.0" rather than exactly 7.19) to accommodate
expected variation from database updates while still catching fundamental errors.
