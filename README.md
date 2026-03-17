# Computational Evidence for De Novo Gene Birth in *Pandoravirus salinus*

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.19046142.svg)](https://doi.org/10.5281/zenodo.19046142)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

This repository contains the complete analytical pipeline and data for:

> **Bradford E.** Computational evidence for a distributed AT-based regulatory mechanism enabling de novo gene birth in *Pandoravirus salinus*. (2026). [Manuscript in preparation]

---

## Summary

*Pandoravirus salinus* has the largest known viral genome (2.47 Mb), with 92.5% of its genes having no detectable homologs outside the Pandoraviridae. This study presents a 12-module computational framework that characterizes *how* new genes emerge from non-coding sequence in this extraordinary genome.

**Key findings:**

- A **five-tier proto-gene continuum** with monotonic gradients in GC content, gene length, and codon usage — from mature genes through evolutionary intermediates to proto-gene candidates
- A **distributed AT-based regulatory system** that pre-loads intergenic regions with promoter and terminator signals (+10.3 / +9.8 percentage points AT enrichment at gene boundaries), resolving how proto-genes can be transcribed regardless of genomic position
- **81.8% of annotated genes** and **97.0% of hypothetical genes** are Pandoravirus-specific orphans, with a 7.2-fold enrichment in external homology for annotated vs. hypothetical genes (OR = 7.19, 95% CI: 3.5–13.6, p = 7.79 × 10⁻¹²)
- A **four-stage maturation ordering**: regulatory signals (pre-existing) → protein structure → codon optimization → GC composition

These findings provide independent computational support for the de novo gene birth hypothesis proposed by the Claverie group (Legendre et al. 2018; Jeudy et al. 2019), arrived at through single-genome compositional analysis rather than multi-strain comparative genomics.

---

## Repository Structure

```
pandoravirus-2.0/
├── README.md                       # This file
├── LICENSE                         # MIT License
├── requirements.txt                # Python dependencies
├── environment.yml                 # Conda environment specification
├── statistics_audit.json           # Maps every manuscript number to its source computation
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
│   └── module12_integration_figures.py
│
├── data/
│   ├── genomes/                    # Source genome files (or download instructions)
│   ├── gff/                        # RNA-seq-informed annotations from IGS GBrowse
│   └── proteins/                   # Exported FASTA files for BLAST
│
├── results/
│   ├── phase1b_results.json        # Transcriptomic validation results
│   ├── phase2_results.json         # BLASTp homology analysis (summary)
│   ├── blastp_annotated.json       # Raw BLAST results (528 annotated proteins)
│   └── blastp_hypothetical.json    # Raw BLAST results (300 hypothetical proteins)
│
├── figures/
│   ├── figure1_continuum.png       # Five-tier proto-gene continuum
│   ├── figure2_regulatory.png      # AT regulatory system + validation
│   ├── figure3_homology.png        # BLASTp homology landscape
│   ├── figure4_maturation.png      # Annotated vs. hypothetical gradient
│   └── supplementary/
│
└── manuscript/
    └── paper_draft.md              # Current manuscript draft
```

---

## Quick Start

### Prerequisites

- Python 3.11+
- PostgreSQL 15+
- Internet access (required for NCBI BLAST queries in Module 11)

### Setup

```bash
# Clone the repository
git clone https://github.com/mathamagician/pandoravirus-2.0.git
cd pandoravirus-2.0

# Create environment
conda env create -f environment.yml
conda activate pandoravirus

# Or use pip
pip install -r requirements.txt
```

### Running the Pipeline

```bash
# Module 1: Download genomes and set up database
python scripts/module01_data_acquisition.py

# Modules 2-10: Compositional analysis and validation (~1-2 hours)
python scripts/module02_genome_architecture.py
python scripts/module03_null_models.py
python scripts/module04_codon_usage.py
python scripts/module05_gene_prediction.py
python scripts/module06_evolutionary_classification.py
python scripts/module07_gene_families.py
python scripts/module08_protein_structure.py
python scripts/module09_regulatory_signals.py
python scripts/module10_transcriptomic_validation.py

# Module 11: BLASTp homology analysis (~24-48 hours, resumable)
python scripts/module11_blastp_homology.py --blast
python scripts/module11_blastp_homology.py --blast-hyp

# Module 12: Generate figures and integrated results
python scripts/module12_integration_figures.py
```

Module 11 (BLASTp) queries NCBI remotely and is fully resumable — it saves results after each query and picks up where it left off if interrupted.

---

## Pipeline Overview

| Module | Description | Key Output |
|--------|-------------|------------|
| 1 | Data acquisition & database setup | PostgreSQL database with 4 genomes |
| 2 | Genome architecture analysis | GC profiles, coding density, intergenic regions |
| 3 | Null model calibration | Markov-1 k-mer expectations, effect-size filtering |
| 4 | Codon usage analysis | Per-gene RSCU, GC3, ENC, CAI |
| 5 | Gene prediction validation | Prodigal vs. GenBank comparison, 347 novel ORFs |
| 6 | Evolutionary classification | GMM clustering, native zone (1.3–2.0 Mb) |
| 7 | Gene families & proto-genes | Five-tier continuum (T1–T5), codon disambiguation |
| 8 | Protein structure prediction | ESM-2 embeddings, structural family clustering |
| 9 | Regulatory signal discovery | Poly-A/T inventory, meta-gene profiles, boundary grammar |
| 10 | Transcriptomic validation | Strand-aware boundary test, hairpin null model |
| 11 | BLASTp homology analysis | 528 annotated + 300 hypothetical proteins vs. NCBI nr |
| 12 | Integration & figures | Manuscript figures 1–4, statistics audit |

---

## Key Results

### Proto-Gene Continuum (Table 1)

| Tier | Description | n | Mean GC | Mean Length |
|------|------------|---|---------|------------|
| T1 | Annotated, Prodigal-confirmed | 516 | 0.650 | 1,452 bp |
| T2 | Hypothetical, Prodigal-confirmed | 784 | 0.648 | 1,082 bp |
| T3 | GenBank-only | 130 | 0.615 | 555 bp |
| T4 | Prodigal high-confidence | 335 | 0.560 | 501 bp |
| T5 | Prodigal low-confidence | 527 | 0.547 | 216 bp |

Spearman ρ = −1.0 for both GC and length across tiers.

### Annotated vs. Hypothetical Homology

| Category | n | Non-self Hit Rate | Orphan Rate |
|----------|---|-------------------|-------------|
| Annotated | 521 | 18.2% | 81.8% |
| Hypothetical | 299 | 3.0% | 97.0% |

Fisher's exact test: OR = 7.19, 95% CI: 3.5–13.6, p = 7.79 × 10⁻¹²

### E-value Robustness

| Threshold | Non-self Hits | Orphan Rate |
|-----------|---------------|-------------|
| 1e-5 | 95 (18.2%) | 81.8% |
| 1e-10 | 88 (16.9%) | 83.1% |
| 1e-20 | 83 (15.9%) | 84.1% |
| 1e-50 | 50 (9.6%) | 90.4% |

---

## Reproducibility Notes

**Deterministic modules (1–10, 12):** Given the same genome sequence and software versions, these modules produce identical results.

**BLAST modules (11):** Results depend on the NCBI nr database, which grows continuously. Exact hit counts may vary with database updates. The key results (orphan rates, OR values, p-values) are robust to small database changes. The random seeds for the hypothetical subsample are fixed (seed 42 for initial 100, seed 43 for additional 200) to ensure the same proteins are selected.

**Statistics audit:** The `statistics_audit.json` file maps every number cited in the manuscript to its source computation and data file, enabling automated verification.

---

## Data Sources

| Data | Source | Accession / URL |
|------|--------|----------------|
| *P. salinus* genome | NCBI GenBank | [NC_022098.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_022098.1) |
| PhiX174 genome | NCBI GenBank | [NC_001422.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_001422.1) |
| Lambda genome | NCBI GenBank | [NC_001416.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_001416.1) |
| Mimivirus genome | NCBI GenBank | [NC_014649.1](https://www.ncbi.nlm.nih.gov/nuccore/NC_014649.1) |
| RNA-seq-informed GFF | IGS GBrowse | [igs.cnrs-mrs.fr/pandoraviruses/](http://www.igs.cnrs-mrs.fr/pandoraviruses/) |
| *A. castellanii* codon table | Kazusa | [Species 5755](https://www.kazusa.or.jp/codon/) |
| BLASTp database | NCBI nr | Accessed March 2026 |

---

## Citation

If you use this pipeline or data, please cite:

```
Bradford E. (2026). Computational evidence for a distributed AT-based regulatory
mechanism enabling de novo gene birth in Pandoravirus salinus. Submitted to Genome Biology and Evolution.
```

This work builds on and provides independent support for:

- Legendre M, et al. (2018). Diversity and evolution of the emerging Pandoraviridae family. *Nature Communications* 9:2285.
- Jeudy S, et al. (2019). *Pandoravirus celtis* illustrates the microevolution processes at work in the giant Pandoraviridae genomes. *Frontiers in Microbiology* 10:430.

---

## Use of AI Tools

AI tools (Anthropic's Claude) were used to assist with code development, data-processing workflow design, analytical planning, and editorial refinement. All analyses, statistical methods, interpretations, and final manuscript content were reviewed and validated by the author.

---

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.
