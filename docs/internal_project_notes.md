# Pandoravirus 2.0 — Transcriptomic Validation, Homology Analysis & Publication

**Analyst:** Eddie Bradford (mathamagician) with Claude Code
**Started:** March 2026
**Builds on:** Pandoravirus 1.0 (10 phases, 5 peer reviews, complete)

## Project Goal

Extend the Pandoravirus 1.0 computational analysis into a publishable result by:
1. Validating the poly-A/T regulatory signal predictions against experimental RNA-seq data
2. Running BLASTp homology analysis on annotated genes to determine evolutionary origins
3. Writing a short research paper on the de novo gene birth model with asymmetric maturation ordering

## Context from Pandoravirus 1.0

**Core thesis:** Pandoravirus salinus is actively generating new genes from non-coding sequence,
enabled by a distributed AT-based regulatory system that pre-loads intergenic regions with
transcriptional capability.

**Key predictions to validate:**
- P1: Comparative genomics — P. dulcis/neocaledonia should show proto-genes at different maturation stages
- **P2: Transcriptomics — poly-A/T runs should coincide with actual TSS/TTS sites** ← Phase 1 target
- P3: AlphaFold — 23 ORFan-domain structural connections should be confirmed
- P4: Experimental — proto-gene candidates should show weak but detectable expression

## Critical Literature Findings

### Already Known (changes our framing)
The Claverie group has ALREADY proposed de novo gene creation from intergenic space as the
main diversifying force of Pandoraviridae (Legendre et al. 2018, Jeudy et al. 2019).
Our contribution is NOT the de novo gene birth model itself, but:
- **Independent confirmation** via compositional/structural analysis (different methodology)
- **The mechanism**: distributed AT regulatory system that pre-loads intergenic space
- **The maturation order**: regulatory (pre-exists) → structural → codon optimization → GC composition
- **The proto-gene continuum**: monotonic gradient across 5 confidence tiers

### Key Papers
1. **Legendre et al. 2018** — "Diversity and evolution of the emerging Pandoraviridae family"
   (Nature Communications 9:2285) — RNA-seq transcriptomics, hairpin termination signals, spliceosomal introns
2. **Jeudy et al. 2019** — "Pandoravirus celtis illustrates the microevolution processes"
   (Frontiers in Microbiology 10:430) — de novo gene creation evidence, open pangenome
3. **Abergel et al. 2023** — "Evolution of giant pandoravirus revealed by CRISPR/Cas9"
   (Nature Communications 14:428) — ORFan genes at 3' end, essential genes at 5' end (P. neocaledonia)
4. **Christo-Foroux et al. 2020** — "A large open pangenome and a small core genome for
   giant pandoraviruses" (Frontiers in Microbiology 11:1486) — open pangenome, ~60 new proteins per strain
5. **Poirot et al. 2019** — "AGCT tetranucleotide anomaly in Pandoravirus genomes"
   (Journal of Virology) — complete absence of AGCT 4-mer, novel DNA editing/selection process
6. **Sun et al. 2015** — "MITE transposable elements in P. salinus genome"
   (BMC Biology) — MITEs colonizing intergenic space, relevant to poly-A/T signal landscape

---

## Phase 1: Transcriptomic Validation — COMPLETE (2026-03-12)

**Goal:** Compare computationally predicted poly-A/T regulatory signals with experimentally
determined transcript boundaries from published RNA-seq data.

**Data limitation:** Raw RNA-seq reads not available in SRA — only accessible via GBrowse
(igs.cnrs-mrs.fr/pandoraviruses/). Analysis used GFF annotations informed by RNA-seq.

### Phase 1a (initial run — had methodological issues)
- Script: `scripts/phase1_transcriptomic_validation.py`
- Hairpin detection overcalibrated at 99.5% (Legendre reported ~70%)
- No 5'/3' ORFan asymmetry in P. salinus (p=0.116) — expected from CRISPR work on P. neocaledonia
- Boundary enrichment not confirmed (wrong test — not strand-aware)

### Phase 1b (corrected — three fixes applied)
- Script: `scripts/phase1b_corrected_analysis.py`
- Figure: `figures/phase1b_corrected.png` (6-panel)
- Results: `results/phase1b_results.json`

**FIX 1 — Intergenic filter:** Removed 2,325 coding runs (16.6%), retained 11,640 intergenic runs.

**FIX 2 — Strand-aware boundary test (HEADLINE RESULT):**
- Promoter signal: 2,332 runs in correct strand-aware position
- Terminator signal: 1,227 runs in correct position
- **74% of boundary-associated runs in correct regulatory position (binomial p < 10⁻⁶)**
- Meta-gene AT spikes: **+10.3pp at starts, +9.8pp at ends** (stronger than Phase 9's +7.0/+9.5pp)

**FIX 3 — Hairpin null model:**
- Observed hairpins: 11.6%, null model: 9.7%, excess only 2.0pp
- Statistically significant (z=3.6) but biologically tiny
- Hairpins NOT the dominant signal — poly-A/T model stands independently
- Legendre's 70% likely used energy-based RNA folding (RNAfold/mfold), not simple palindrome matching

**Two-population test (active vs latent regulatory substrate):**
- Boundary-proximal: 4,331 runs; Dispersed: 7,309 runs
- NO compositional difference (A/T ratio p=0.105, length p=0.093)
- Supports "latent regulatory substrate" model — dispersed runs are compositionally identical to boundary-active ones

### Phase 1 Conclusions
1. Strand-aware poly-A/T regulatory signal **VALIDATED and strengthened**
2. Hairpin termination (Legendre 2018) operates at a different analytical level — does not contradict poly-A/T model
3. No 5'/3' ORFan asymmetry in P. salinus (uniform distribution supports "gene birth everywhere")
4. Dispersed poly-A/T runs = latent regulatory substrate (same composition as active boundary runs)
5. Prediction 2 partially validated (full validation requires raw RNA-seq TSS/TTS coordinates)

---

## Phase 2: BLASTp Homology Analysis — COMPLETE (2026-03-14)

**Goal:** Determine evolutionary origins and taxonomic affinities of the 528 annotated
(non-hypothetical) P. salinus genes via BLASTp against NCBI nr.

- Script: `scripts/phase2_blastp_homology.py`
- Figure: `figures/phase2_blastp_homology.png` (5-panel)
- Results: `results/phase2_results.json` (analysis summary), `results/blastp_annotated.json` (raw BLAST, 3.5 MB)
- Data: `data/proteins_annotated.fasta` (528), `data/proteins_hypothetical.fasta` (902)
- Parameters: E-value < 1e-5, top 10 hits per protein, searched against nr
- Runtime: ~2 days via remote NCBI BLAST (fully resumable)

### Headline: 81.8% of annotated genes have ONLY Pandoravirus hits

| Category | Count | % of successful |
|---|---|---|
| Total attempted | 528 | — |
| Successful BLAST queries | 521 | 100% denominator |
| Connection errors (excluded) | 7 | — |
| With any BLAST hit | 521 | 100.0% |
| **Pandoravirus-only hits** | **426** | **81.8%** |
| With non-self (non-Pandoravirus) hit | 95 | 18.2% |
| Genuine no hits | 0 | 0.0% |

*Note: 7 "no hits" in original run were connection failures (timeouts), not genuine zero-hit results.*

### Taxonomy of non-self top hits (n=95) — CORRECTED

*Original analysis misclassified 12 marine green algae as "Bacteria" due to substring match
on 'monas' in organism names. Fixed 2026-03-14 with curated organism-level lookup table.*

| Category | Count | % of non-self | Mean % Identity |
|---|---|---|---|
| Marine algae (prasinophytes + brown algae) | 28 | 29.5% | 44.2% |
| Giant virus (NCLDVs) | 24 | 25.3% | 47.5% |
| Amoebozoa (host lineage) | 13 | 13.7% | 50.1% |
| Bacteria (genuine prokaryotes) | 11 | 11.6% | 38.1% |
| Insect (aphids — likely Wolbachia-mediated ankyrin) | 11 | 11.6% | 41.1% |
| Other eukaryote (fungi, oomycetes, rotifers, etc.) | 7 | 7.4% | 42.9% |
| Other virus (baculovirus) | 1 | 1.1% | — |

### Identity distribution (non-self top hits)
- Mean: 44.7%, Median: 42.3%, Min: 25.8%, Max: 84.0%
- <30% identity: 5 (5.3%) — very distant
- 30-50% identity: 65 (68.4%) — distant homology
- >50% identity: 25 (26.3%) — moderate-to-clear homology

### Spatial analysis — native zone enrichment
- Non-self hits **enriched in 1.3-2.0 Mb native zone**: OR=1.95, **Fisher's p=0.004**
- 45 of 95 non-self hits fall in native zone (expected ~27 by genome proportion)
- Confirms native zone contains evolutionarily older genes with broader phylogenetic connections

### Product-level hit rates
| Product | Count | Non-self hit rate | Mean identity | Top taxonomy |
|---|---|---|---|---|
| Ankyrin repeat | 130 | 42.3% | 44.2% | Marine algae |
| F-box domain | 92 | 8.7% | 42.1% | Bacteria |
| MORN repeat | 44 | 4.5% | 39.8% | Bacteria |
| F-box incomplete | 26 | 7.7% | 45.4% | Amoebozoa |
| Ubiquitin | 5 | 40.0% | 46.6% | Other virus |
| Fascin-like, Atrophin, DHFR, Ring | 29 | 0% | — | Pandoravirus-specific |

*The product-level gradient — from "clearly acquired" (ankyrin 42.3%) through "diverged beyond
recognition" (F-box 8.7%) to "functionally convergent or de novo" (Fascin 0%) — is direct
evidence for different evolutionary ages of domain acquisition.*

### Phase 2 Conclusions
1. **81.8% of annotated genes are unique to Pandoraviruses** — strong de novo origin support
2. Non-self homologs show **eukaryote-dominated** pattern: marine algae (29.5%) + amoebozoa (13.7%) = 43%. Shared aquatic ecology likely drives gene exchange — not generic multi-kingdom HGT
3. Native zone enrichment (p=0.004) connects to Phase 6 maturation ordering — older genome region retains more external homologs
4. Most homology is distant (median 42.3%) — ancient divergence or convergent domain evolution
5. Product-level gradient from "clearly acquired" to "Pandoravirus-specific" is the strongest evidence for different evolutionary ages of domain acquisition

---

## Phase 2b: Hypothetical Gene Subsample — COMPLETE (2026-03-15)

**Goal:** Compare non-self hit rates between annotated and hypothetical genes to establish
the maturation gradient predicted by the de novo gene birth model.

- Subsample: 300 random hypothetical proteins (seed=42 for first 100, seed=43 for additional 200)
- 299 successful queries, 1 connection error, 1 genuine no-hit

### Headline: 97.0% of hypothetical genes are Pandoravirus-only

| Category | Annotated (n=521) | Hypothetical (n=299) |
|---|---|---|
| **Non-self hit rate** | **18.2%** | **3.0%** |
| Pandoravirus-only | 81.8% | 97.0% |
| Non-self hits | 95 | 9 |
| Fisher's exact | — | **OR=7.19, 95% CI: 3.5–13.6, p=7.79e-12** |

### Hypothetical non-self hit taxonomy (n=9)

| Category | Count | Interpretation |
|---|---|---|
| Giant virus (NCLDVs) | 6 (66.7%) | Shared NCLDV gene pool |
| Amoebozoa (host) | 2 (22.2%) | Recent host→virus HGT |
| Bacteria | 1 (11.1%) | Distant homology |
| Marine algae | 0 | — |
| Insect | 0 | — |

*The hypothetical hits connect exclusively to phylogenetically proximate sources (NCLDV pool +
host lineage). The more distant connections (marine algae, insects) appear only in annotated
genes — a second gradient supporting the maturation model.*

- Hypothetical spatial enrichment: 7/9 hits in native zone (OR=11.9, p=0.0009)

### Three-level maturation gradient

| Gene Category | n | Non-self Hit Rate | Interpretation |
|---|---|---|---|
| Annotated (with recognized domains) | 521 | 18.2% | Oldest genes; some retain detectable homology |
| Hypothetical (ORFs, no domain annotation) | 299 | 3.0% | Younger genes; almost entirely orphans |
| Intergenic regions | ~1,900 IGRs | 0% (by definition) | Not yet genes; the raw substrate |

*This gradient directly supports the de novo gene birth model with maturation staging.*

### E-value Sensitivity Analysis

| Threshold | Non-self hits | Hit rate | Orphan rate |
|---|---|---|---|
| < 1e-5 (baseline) | 95 | 18.2% | 81.8% |
| < 1e-10 | 88 | 16.9% | 83.1% |
| < 1e-20 | 83 | 15.9% | 84.1% |
| < 1e-50 | 50 | 9.6% | 90.4% |

*Result robust across thresholds — graceful degradation confirms signal is not threshold-dependent.*

---

## Phase 3: Paper Writing — NOT STARTED

**Goal:** Write a short research paper for a computational biology or virology journal.

**Target framing:** Independent computational confirmation of de novo gene birth in Pandoravirus,
with novel discovery of asymmetric maturation ordering and distributed AT regulatory mechanism.

**Three novel contributions:**
1. Maturation ordering: regulatory signals pre-exist → structural features evolve → codon optimization → GC composition
2. AT regulatory mechanism: distributed poly-A/T system that pre-loads intergenic space with regulatory capability
3. Comparative-free pipeline: all findings derived from single-genome analysis without requiring multiple strain comparison

**Paper structure:**
- Introduction: Pandoravirus genome, ORFan problem, existing de novo gene birth hypothesis
- Methods: Pipeline description (10 phases, comparison ladder, null models)
- Results: Proto-gene continuum, maturation order, regulatory signal discovery, transcriptomic validation, BLASTp homology
- Discussion: Relationship to Claverie group findings, mechanism contribution, predictions
- Supplementary: Full pipeline code, all figures, database schema

**Peer review framing guidance:** See `feedback-pandoravirus2-review.md` in Claude Code memory.

---

## File Inventory

### Scripts
| File | Description | Status |
|---|---|---|
| `scripts/phase1_transcriptomic_validation.py` | Phase 1a (initial, has issues) | Superseded |
| `scripts/phase1b_corrected_analysis.py` | Phase 1b (corrected, three fixes) | Complete |
| `scripts/phase2_blastp_homology.py` | Phase 2 BLASTp analysis | Complete |

### Figures
| File | Description |
|---|---|
| `figures/phase1_transcriptomic_validation.png` | Phase 1a (4-panel, superseded) |
| `figures/phase1b_corrected.png` | Phase 1b (6-panel, current) |
| `figures/phase2_blastp_homology.png` | Phase 2 (7-panel: A-E + F annotated-vs-hyp + G maturation gradient) |

### Results
| File | Description | Size |
|---|---|---|
| `results/phase1_results.json` | Phase 1a raw results | 2 KB |
| `results/phase1b_results.json` | Phase 1b corrected results | 1 KB |
| `results/blastp_annotated.json` | Phase 2 raw BLAST (all 528 proteins) | 3.5 MB |
| `results/blastp_hypothetical.json` | Phase 2b raw BLAST (100 hypothetical proteins) | ~500 KB |
| `results/phase2_results.json` | Phase 2 analysis summary (incl. comparison + CIs) | ~115 KB |

### Data
| File | Description |
|---|---|
| `data/proteins_annotated.fasta` | 528 annotated protein sequences |
| `data/proteins_hypothetical.fasta` | 902 hypothetical protein sequences |

### Docs
| File | Description |
|---|---|
| `docs/literature_review.md` | Full literature review |

## Tech Stack
- Python 3.11.9, BioPython 1.86, psycopg2, numpy, scipy, scikit-learn, matplotlib
- PostgreSQL 15 (pandoravirus database, localhost:5432)
- NCBI BLAST (remote via BioPython NCBIWWW)
- Existing Phase 9 regulatory signal data as baseline

## Original Project
- Location: `C:\Users\Eddie\Code Projects\Pandoravirus`
- Database: pandoravirus (PostgreSQL, localhost:5432)
- All 10 phases complete with 35 key findings
