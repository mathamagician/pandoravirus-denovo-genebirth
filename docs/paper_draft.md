# Computational evidence for a distributed AT-based regulatory mechanism enabling de novo gene birth in *Pandoravirus salinus*

**Eddie Bradford**

---

## Abstract

*Pandoravirus salinus* possesses the largest known viral genome (2.47 Mb) with an unprecedented 92.5% of genes classified as ORFans — proteins with no detectable homologs outside the Pandoraviridae. While de novo gene creation from non-coding sequence has been proposed as the primary diversifying mechanism of this lineage, the genomic infrastructure enabling this process remains poorly understood. Here we present a computational framework combining compositional analysis, structural prediction, transcriptomic validation, and genome-wide homology analysis to characterize the mechanism underlying de novo gene birth in *P. salinus*.

We identify a five-tier proto-gene continuum exhibiting a monotonic gradient from established genes (mean GC 0.650, 1,452 bp) through proto-genes approaching intergenic composition (GC 0.547, 216 bp), with codon optimization detectable even in the most immature candidates (p = 0.0005). Critically, we discover that poly-A/T regulatory signals — functioning as promoters and terminators — pre-exist throughout intergenic space, with proto-genes exhibiting mature regulatory signatures indistinguishable from established genes (p = 0.14). This distributed AT-based regulatory system produces the strongest gene boundary signals observed across all genomes tested (+10.3 and +9.8 percentage points AT enrichment at gene starts and ends, respectively), despite operating against the genome's 62% GC composition.

Strand-aware analysis of poly-A/T positioning relative to experimentally informed gene annotations confirms that 74% of boundary-associated runs occupy the correct regulatory position (binomial p < 10⁻⁶). BLASTp analysis reveals that 81.8% of annotated genes and 97.0% of hypothetical genes are Pandoravirus-specific orphans, with annotated genes 7.2-fold more likely to retain external homologs than hypothetical genes (Fisher's exact test, OR = 7.19, 95% CI: 3.5–13.6, p = 7.79 × 10⁻¹²). This gradient is robust across E-value thresholds (orphan rate 81.8%–90.4% from 10⁻⁵ to 10⁻⁵⁰).

These findings establish a four-stage maturation model — regulatory signals (pre-existing) → protein structure → codon optimization → GC composition — and provide an independent computational line of support, via different methodology, for the de novo gene birth hypothesis proposed by the Claverie group. The distributed AT regulatory system resolves a key puzzle: proto-genes need not arise near existing promoters because transcriptional infrastructure pervades the intergenic landscape.

---

## Introduction

The discovery of *Pandoravirus salinus* revealed the largest viral genome known to science — 2.47 megabases encoding approximately 1,430 predicted proteins, of which 92.5% have no detectable homologs outside the Pandoraviridae (Philippe et al. 2013). This unprecedented ORFan fraction poses a fundamental evolutionary question: where do these genes come from?

The Claverie group has provided compelling evidence that de novo gene creation from non-coding sequence is the primary diversifying mechanism of the Pandoraviridae. Legendre et al. (2018) conducted strand-specific RNA-seq across four Pandoravirus strains, demonstrating that 83–87% of genome sequence is transcribed, that palindromic hairpin structures mark approximately 70% of transcript 3' ends, and that the Pandoravirus pangenome is open — each new strain adds roughly 60 previously unseen proteins. Jeudy et al. (2019) provided direct evidence by comparing the nearly identical genomes of *P. celtis* and *P. quercus* (96.7% identity), identifying genes in one strain that correspond to intergenic regions in the other — genes captured in the act of being born. Aherfi et al. (2018) confirmed the open pangenome with a small conserved core, and Bisio et al. (2023) demonstrated through CRISPR/Cas9 knockout that essential genes concentrate at the 5' end while ORFans accumulate at the 3' end of *P. neocaledonia*.

These studies establish *what* is happening — genes are emerging from non-coding sequence — but leave a mechanistic gap: *how* does intergenic space acquire the regulatory and coding capacity to produce functional genes? If new genes require promoters and terminators to be expressed, how do they arise in intergenic regions distant from existing regulatory elements? And if gene birth is a gradual process, what is the sequence of events from non-coding DNA to functional protein?

Additional compositional analyses have revealed that Pandoravirus genomes are subject to unusual sequence constraints. Poirot et al. (2019) documented the near-complete absence of the AGCT tetranucleotide in some strains, pointing to a novel DNA editing or selection process. Sun et al. (2015) identified miniature inverted-repeat transposable elements (MITEs) colonizing the intergenic space of *P. salinus*, demonstrating that multiple forces shape the non-coding landscape from which new genes emerge.

Here we present a multi-method computational framework that addresses the mechanistic gap. Through a 12-module analytical pipeline encompassing genome architecture, null model calibration, codon usage analysis, gene prediction validation, evolutionary classification, gene family detection, protein structure prediction, regulatory signal discovery, transcriptomic validation, and genome-wide homology analysis, we provide three novel contributions. First, we identify a quantitative five-tier proto-gene continuum with monotonic gradients across all measured features, providing a maturation metric for gene birth. Second, we discover a distributed AT-based regulatory system that pre-loads intergenic regions with promoter and terminator signals, resolving the spatial puzzle of how genes arise far from existing regulatory elements. Third, we establish an asymmetric maturation ordering — regulatory signals pre-exist, followed by structural features, codon optimization, and finally GC composition — that defines the temporal sequence of de novo gene birth.

Our approach differs fundamentally from previous work: all findings derive from single-genome computational analysis of *P. salinus* without requiring multi-strain comparison, providing an independent line of evidence that complements the comparative genomics approach of the Claverie group.

---

## Methods

### Genome and database

The complete genome of *Pandoravirus salinus* (GenBank accession NC_022098.1; 2,473,870 bp, 1,430 annotated CDS) was loaded into a PostgreSQL 15 database alongside three comparison genomes spanning the dsDNA virus size spectrum: bacteriophage PhiX174 (NC_001422.1; 5,386 bp), bacteriophage Lambda (NC_001416.1; 48,502 bp), and *Acanthamoeba polyphaga* mimivirus (NC_014649.1; 1,181,549 bp). The database schema comprises 22 tables including position-indexed nucleotide sequences, gene annotations, intergenic regions, codon usage statistics, and derived analytical tables. All analyses used Python 3.11 with BioPython 1.86, NumPy, SciPy, scikit-learn, and matplotlib.

### Compositional analysis pipeline (Modules 2–7)

A 12-module analytical pipeline was applied to characterize the compositional and structural landscape of *P. salinus*. Here we summarize the methods underlying findings presented in Results; full pipeline details and scripts are available in the Supplementary Materials.

**Genome architecture (Module 2).** Eight metrics were computed genome-wide and per-gene: GC content, coding density, strand bias, intergenic region length, Shannon entropy, gene length distribution, and sliding-window GC profiles. All metrics were compared across the four-genome ladder.

**Null model calibration (Module 3).** Markov-1 null models preserving dinucleotide frequencies were used to assess k-mer enrichment and depletion (k = 4, 6, 8). Effect-size filtering (observed/expected < 0.5 or > 2.0) was applied alongside z-score significance to reduce false positives in large genomes, where 88.6% of z > 3 signals lacked meaningful effect sizes.

**Codon usage analysis (Module 4).** Per-gene relative synonymous codon usage (RSCU), GC3 (wobble position), GC12 (constrained positions), effective number of codons (ENC), and codon adaptation index (CAI) were computed. CAI was calculated using *Acanthamoeba castellanii* codon frequencies from the Kazusa database. Annotated versus ORFan gene comparisons tested whether ORFans show distinct compositional signatures indicative of foreign acquisition.

**Gene prediction validation (Module 5).** Ab initio gene prediction with Prodigal was compared against GenBank annotations using ≥80% reciprocal overlap matching. Novel Prodigal-only predictions (n = 347) were characterized by length, GC content, and genomic position.

**Evolutionary origin classification (Module 6).** Gaussian Mixture Model (GMM) clustering on five per-gene features (GC, GC3, GC12, ENC, CAI) identified compositional clusters. The native zone (1.3–2.0 Mb), defined independently by this clustering, was used in subsequent spatial analyses.

**Gene family and proto-gene analysis (Modules 7–7c).** All-vs-all protein 3-mer frequency cosine similarity (1,021,735 pairwise comparisons) with hierarchical clustering identified gene families. Null calibration using composition-matched random proteins confirmed that the diffuse similarity cloud (mean 0.109) is a compositional artifact (null mean 0.105; Cohen's d = 0.06), while the 14 detected families represent genuine signals.

Proto-gene candidates were classified into five confidence tiers (T1–T5) based on prediction agreement between GenBank annotations and Prodigal: T1, annotated genes confirmed by both (n = 516); T2, hypothetical genes confirmed by both (n = 784); T3, GenBank-only genes missed by Prodigal (n = 130); T4, high-confidence Prodigal-only predictions (n = 335); T5, low-confidence Prodigal-only predictions (n = 527). Disambiguation tests compared T4/T5 codon usage against length-matched random intergenic ORFs and assessed spatial proximity to established genes.

**Protein structure prediction (Module 8).** ESM-2 (650M parameters; Lin et al. 2023) embeddings were generated for 170 representative proteins spanning all tiers and compositional clusters. PCA reduction (top 20 components, 85.8% variance) followed by hierarchical clustering identified structural families. Null calibration confirmed structural signal genuineness (within-family similarity 0.533 vs null 0.368; Cohen's d = 1.45), though 41% of PC1 variance correlates with GC content.

### Regulatory signal analysis (Module 9)

Seven complementary analyses characterized the regulatory landscape: (i) genome-wide poly-A/T run inventory (runs ≥ 5 bp), (ii) meta-gene AT content profiles across gene boundaries, (iii) upstream k-mer enrichment analysis, (iv) intergenic grammar comparison between convergent and divergent regions, (v) strand-asymmetric boundary positioning, (vi) ORFan versus annotated regulatory signature comparison, and (vii) proto-gene regulatory context assessment.

### Transcriptomic validation (Module 10)

Computationally predicted regulatory signals were compared against experimentally informed gene annotations from Legendre et al. (2018). GFF annotations incorporating RNA-seq evidence were used in lieu of raw reads, which are not deposited in public repositories.

Three methodological corrections were applied relative to an initial analysis: (i) an intergenic filter removed 2,325 poly-A/T runs (16.6%) overlapping coding sequences, retaining 11,640 intergenic runs; (ii) a strand-aware boundary test classified each run as promoter-positioned (poly-A upstream of sense-strand gene start, or poly-T upstream of antisense-strand gene start) or terminator-positioned, replacing a strand-naive proximity test; (iii) a hairpin null model using strict parameters (stem ≥ 6 bp, exact complementarity) quantified the contribution of palindromic sequences to termination signals.

The two-population test compared AT composition and length distributions between boundary-proximal runs (within 200 bp of a gene boundary; n = 4,331) and dispersed intergenic runs (n = 7,309) to test whether boundary-active and dispersed runs represent distinct populations or a uniform regulatory substrate.

### BLASTp homology analysis (Modules 11–12)

All 528 annotated (non-hypothetical) *P. salinus* protein sequences were searched against the NCBI non-redundant protein database (nr) using remote BLASTp via BioPython (E-value threshold 1e-5, top 10 hits per query). Of 528 queries, 521 returned valid results; 7 connection failures were excluded from analysis. A subsample of 300 randomly selected hypothetical proteins (from 902 total; random seed 42 for initial 100, seed 43 for additional 200) was searched under identical parameters, yielding 299 successful queries.

Non-self hits were defined as hits to organisms other than Pandoravirus strains. Taxonomic classification used a curated organism-level lookup table covering all 50 unique organisms encountered in the dataset, supplemented by an ordered keyword fallback and NCBI Entrez API lookup for unrecognized organisms. During quality auditing, the initial automated classification was found to have misidentified 12 prasinophyte green algae as bacteria due to a substring match on the suffix '-monas'; this was corrected by implementing the curated lookup table approach.

Identity distribution, spatial enrichment in the native zone (1.3–2.0 Mb), and product-level hit rates were computed. The annotated versus hypothetical comparison used Fisher's exact test with Haldane-corrected Woolf method for the 95% confidence interval on the odds ratio, and Wilson score intervals for proportions. E-value sensitivity analysis re-filtered existing hits at thresholds of 10⁻¹⁰, 10⁻²⁰, and 10⁻⁵⁰ to assess robustness.

---

## Results

### A five-tier proto-gene continuum

Classification of all 2,292 predicted ORFs (1,430 GenBank + 862 additional Prodigal predictions) into five confidence tiers revealed a monotonic gradient across every measured compositional feature (Table 1). GC content declined continuously from 0.650 in mature annotated genes (T1) to 0.547 in the most immature proto-gene candidates (T5) — approaching the intergenic average of 0.539. Gene length followed the same trajectory (1,452 bp → 216 bp), as did GC3 (0.750 → 0.672 for the tiers where wobble position could be computed) and CAI (0.792 → 0.762).

**Table 1. Five-tier proto-gene continuum.** Spearman ρ = −1.0 for both GC and length across tiers.

| Tier | Description | n | Mean length (bp) | Mean GC | GC3 | CAI |
|------|------------|---|------------------|---------|-----|-----|
| T1 | Annotated, Prodigal-confirmed | 516 | 1,452 | 0.650 | 0.750 | 0.792 |
| T2 | Hypothetical, Prodigal-confirmed | 784 | 1,082 | 0.648 | 0.745 | 0.787 |
| T3 | GenBank-only (Prodigal missed) | 130 | 555 | 0.615 | 0.672 | 0.762 |
| T4 | Prodigal high-confidence only | 335 | 501 | 0.560 | — | — |
| T5 | Prodigal low-confidence only | 527 | 216 | 0.547 | — | — |

The gradient did not merely reflect annotation quality. Codon usage analysis comparing T4 and T5 predictions against length-matched random intergenic ORFs revealed statistically significant genome-directed optimization even in the least mature candidates: T4 codon usage was significantly closer to the genome-preferred pattern than random (p = 0.020), with T5 showing an even stronger signal (p = 0.0005). The effect size was modest — approximately 10–15% of the distance from random to fully optimized — consistent with weak but detectable selection on recently emerged ORFs. Prodigal does not incorporate codon preference in its prediction algorithm, so this signal must originate from the sequences themselves.

Proto-genes showed no spatial clustering near established genes (T4 p = 0.99; T5 p = 1.0), ruling out regulatory co-option or operon extension as mechanisms. Instead, proto-gene candidates were distributed uniformly across the genome, consistent with a model in which gene birth can occur anywhere in the intergenic landscape.

Compositional uniformity analysis provided additional context. Despite being 460-fold larger than PhiX174, *P. salinus* exhibits the lowest GC coefficient of variation (6.1%) of all tested genomes, and ORFan genes are compositionally indistinguishable from annotated genes across all metrics (GC3 difference +0.015, CAI difference +0.008, upstream 6-mer cosine similarity 0.9824). This homogeneity argues against recent patchwork horizontal gene transfer, which would produce detectable compositional heterogeneity.

### A distributed AT-based regulatory system

The most striking finding of the compositional analysis was the discovery of an AT-based regulatory system operating against the genome's high GC content. After filtering to intergenic runs only (see Transcriptomic validation below), *P. salinus* exhibited gene boundary AT signals dramatically stronger than any other genome tested: +10.3 percentage points AT enrichment at gene starts and +9.8 percentage points at gene ends, compared to approximately +1 percentage point for Lambda and +1–2 percentage points for Mimivirus (Table 2). This signal runs counter to the genome's 62% GC composition, requiring active regulatory selection to maintain.

**Table 2. Gene boundary AT signals across genomes.** *P. salinus* values are from filtered intergenic analysis; unfiltered values (+7.0/+9.5 pp) are lower because coding-region poly-A/T runs dilute the signal.

| Genome | AT content (GC content) | Start ΔAT (pp) | End ΔAT (pp) | Poly-A/T intergenic (%) |
|--------|------------------------|----------------|-------------|------------------------|
| Lambda | 50.1% (49.9% GC) | −1.1 | +1.6 | 12–24% |
| Mimivirus | 72.0% (28.0% GC) | +1.1 | +2.1 | 16–17% |
| **P. salinus** | **38.3% (61.7% GC)** | **+10.3** | **+9.8** | **82.8–83.9%** |

Poly-A and poly-T homopolymer runs (≥5 bp) numbered 6,955 and 7,010 respectively, with 82.8–83.9% located in intergenic regions — a pattern unique to *P. salinus* (Lambda: 12–24%; Mimivirus: 16–17%). These runs displayed strand-asymmetric positioning consistent with regulatory function: poly-A runs were enriched upstream of gene starts (41.5% within 200 bp), while poly-T runs were enriched downstream of gene ends (40.8%).

Divergent intergenic regions (flanked by genes transcribed outward in both directions) were 17% longer than convergent regions (582 bp vs 498 bp; p < 0.0001), consistent with a requirement for two promoters to be accommodated. However, divergent regions were less dense in poly-A runs (6.67 vs 7.69 runs/kb; p = 0.019), indicating that the extra length represents spacer DNA between two normal-density promoter regions.

The critical finding for the de novo gene birth model was that proto-gene candidates exhibited regulatory signals indistinguishable from established genes. T4 upstream AT content was 45.0% (p = 0.14 vs annotated 45.9%) and T5 was 46.7% (p = 0.25) — neither significantly different from mature genes. This means that the regulatory infrastructure for transcription pre-exists in intergenic space before ORFs emerge, resolving the paradox of how spatially random proto-genes can be expressed: they are born into a landscape already equipped with promoters and terminators.

### Transcriptomic validation

To validate the computationally predicted regulatory signals, we compared poly-A/T run positions against experimentally informed gene annotations from Legendre et al. (2018). After filtering to 11,640 intergenic runs and applying strand-aware classification, we found that 74% of boundary-associated runs occupied the correct regulatory position — poly-A upstream on the sense strand or poly-T upstream on the antisense strand for promoters, and the converse for terminators (binomial p < 10⁻⁶).

This strand-aware analysis identified 2,332 runs in promoter-consistent positions and 1,227 in terminator-consistent positions. The asymmetry (1.9:1 promoter-to-terminator ratio) may reflect the greater complexity of transcription initiation relative to termination, or may indicate that a subset of termination events uses the hairpin mechanism described by Legendre et al. (2018) rather than poly-T signals.

Meta-gene AT content profiles — averaging AT content across aligned gene boundaries — revealed sharp spikes at transcript boundaries: +10.3 percentage points at gene starts and +9.8 percentage points at gene ends. These values exceed the initial unfiltered estimates (+7.0/+9.5 pp from Module 9) because the intergenic filter removed coding-region runs that diluted the signal.

To assess whether palindromic hairpin structures, proposed by Legendre et al. (2018) as the primary termination signal, could account for the poly-A/T signal, we implemented a hairpin null model with strict parameters (stem ≥ 6 bp, exact complementarity). The observed hairpin fraction (11.6%) exceeded the null model expectation (9.7%) by only 2.0 percentage points (z = 3.6, statistically significant but biologically small). The 70% hairpin frequency reported by Legendre et al. (2018) likely reflects energy-based RNA folding analysis (e.g., RNAfold/mfold) rather than simple palindrome detection. Importantly, the poly-A/T signal stands independently of hairpin structures.

The two-population test comparing boundary-proximal (n = 4,331) and dispersed intergenic (n = 7,309) poly-A/T runs found no compositional difference in AT ratio (p = 0.105) or length (p = 0.093). This supports the latent regulatory substrate model: dispersed runs throughout intergenic space are compositionally identical to those actively functioning at gene boundaries. The regulatory potential is distributed, not concentrated — any region of the genome could, in principle, support transcription of a newly emerged ORF.

### Homology landscape and the maturation gradient

BLASTp analysis of all 528 annotated *P. salinus* proteins against the NCBI nr database revealed that 426 of 521 successful queries (81.8%) returned hits exclusively to other Pandoravirus strains, confirming the extreme lineage-specificity of this proteome. The 95 proteins with non-self hits showed an eukaryote-dominated taxonomic distribution (Table 3): marine algae (28, 29.5%), giant viruses (24, 25.3%), Amoebozoa (13, 13.7%), bacteria (11, 11.6%), insects (11, 11.6%), other eukaryotes (7, 7.4%), and other viruses (1, 1.1%).

**Table 3. Taxonomy and identity of non-self BLAST hits (annotated genes, n = 95).**

| Category | Count | % | Mean identity |
|----------|-------|---|--------------|
| Marine algae (prasinophytes + brown algae) | 28 | 29.5 | 44.2% |
| Giant virus (NCLDVs) | 24 | 25.3 | 47.5% |
| Amoebozoa (host lineage) | 13 | 13.7 | 50.1% |
| Bacteria | 11 | 11.6 | 38.1% |
| Insect (aphids, likely endosymbiont-mediated) | 11 | 11.6 | 41.1% |
| Other eukaryote | 7 | 7.4 | 42.9% |
| Other virus | 1 | 1.1 | — |

The marine algae category represents an ecological grouping, not a phylogenetic one: green algae (Chlorophyta/Prasinophyceae; *Micromonas*, *Bathycoccus*, *Ostreococcus*) and brown algae (Phaeophyceae; *Ectocarpus*, *Scytosiphon*) are as distantly related as animals and fungi, yet both share marine coastal habitat with *P. salinus*. The 11 insect hits mapped exclusively to ankyrin and F-box domains in aphid genomes (*Acyrthosiphon pisum*, *Aphis* spp.); these likely reflect convergent domain expansion mediated by Wolbachia endosymbionts rather than direct Pandoravirus-to-insect horizontal transfer.

The taxonomic composition of non-self hits identifies the ecological neighborhoods in which Pandoravirus proteins retain detectable connections, but taxonomy alone does not test the maturation model. Non-self hits were spatially enriched in the native zone (1.3–2.0 Mb): 45 of 95 hits fell within this region (OR = 1.94, Fisher's exact p = 0.004), consistent with the native zone containing evolutionarily older genes that retain broader phylogenetic connections. Product-level analysis revealed a gradient from clearly acquired domains (ankyrin repeat: 42.3% non-self hit rate) through diverged domains (F-box: 8.7%) to Pandoravirus-specific domains with zero external homologs (Fascin, Atrophin, DHFR, Ring: 0%), providing direct evidence for different evolutionary ages of domain acquisition.

The more informative question is whether homology tracks gene maturity — the central prediction of the maturation model. We therefore compared annotated and hypothetical proteins directly. The annotated-versus-hypothetical comparison provided the strongest test of the maturation model. Of 299 successfully queried hypothetical proteins, only 9 (3.0%) had non-self hits, compared to 95 of 521 (18.2%) for annotated proteins — a 7.2-fold enrichment (Fisher's exact test, OR = 7.19, 95% CI: 3.5–13.6, p = 7.79 × 10⁻¹²; Table 4). This result was robust across E-value thresholds: the annotated orphan rate remained above 80% even at the most stringent threshold tested (81.8% at 10⁻⁵, 83.1% at 10⁻¹⁰, 84.1% at 10⁻²⁰, 90.4% at 10⁻⁵⁰).

**Table 4. Annotated vs hypothetical homology comparison.**

| Gene category | n | Non-self hit rate | Interpretation |
|--------------|---|------------------|----------------|
| Annotated (recognized domains) | 521 | 18.2% | Oldest genes; some retain detectable homology |
| Hypothetical (ORFs, no domains) | 299 | 3.0% | Younger genes; almost entirely orphans |
| Intergenic regions | ~1,900 IGRs | 0% (by definition) | Not yet genes; the raw substrate |

This comparison reveals not only a quantitative difference in orphan rate, but also a qualitative contraction in phylogenetic reach among younger genes. The 9 hypothetical non-self hits displayed a taxonomic profile distinct from annotated hits: 6 giant viruses (66.7%), 2 Amoebozoa (22.2%), and 1 bacterium (11.1%) — with zero marine algae and zero insects. Hypothetical genes with detectable external homology connect exclusively to phylogenetically proximate sources (the NCLDV gene pool and the amoebozoan host lineage), while more distant connections (marine algae, insects, diverse bacteria) appear only in annotated genes. This phylogenetic distance gradient — younger genes retain similarity only to nearest relatives — is an additional prediction of the maturation model confirmed by the data.

---

## Discussion

### A mechanistic model for de novo gene birth

The results presented here converge on a coherent model for de novo gene birth in *P. salinus* that can be articulated as a three-step logical chain. First, poly-A/T runs with regulatory positioning exist throughout the intergenic space, as demonstrated by the 82.8–83.9% intergenic fraction, the strand-asymmetric boundary enrichment, and the +10.3/+9.8 pp meta-gene AT spikes confirmed by strand-aware validation. Second, boundary-proximal and dispersed intergenic runs are compositionally identical (AT ratio p = 0.105, length p = 0.093), indicating that the regulatory substrate is distributed uniformly rather than concentrated at existing gene boundaries. Third, proto-genes emerge into these pre-existing regulatory contexts with mature transcriptional signals already in place (T4 p = 0.14, T5 p = 0.25), acquire coding features through a gradual maturation process measurable as the five-tier continuum, and progressively accumulate external homologs as they age (3.0% → 18.2% non-self hit rate).

This model predicts a specific maturation ordering. Regulatory signals pre-exist as a feature of intergenic sequence composition, representing the earliest (indeed, prerequisite) stage. Protein structural features emerge next — proto-genes are not structurally isolated but share recognizable folds with established domains (nearest-neighbor structural similarity p = 0.90 vs established genes). Codon optimization follows, detectable as a 10–15% shift from random toward genome-preferred usage (T4 p = 0.020, T5 p = 0.0005). GC composition shifts last, as the slowest and most passive process, with proto-genes (GC 0.547) only gradually approaching the coding average (GC 0.650) from the intergenic baseline (GC 0.539).

An alternative interpretation is that the T4 and T5 proto-gene tiers reflect annotation artifacts or false-positive Prodigal predictions rather than genuine evolutionary intermediates. Several observations argue against this. First, both T4 and T5 exhibit codon usage significantly shifted toward genome-preferred patterns relative to length-matched random intergenic ORFs, with the strongest signal retained even in T5 (p = 0.0005), indicating sequence-level structure inconsistent with random noise. Second, representative proto-genes are not structurally isolated but cluster near established domain families in embedding space, suggesting continuity in protein-level organization rather than arbitrary ORF detection. Third, proto-genes show no spatial clustering near established genes, arguing against simple operon extension or boundary spillover from neighboring annotations. Taken together, these results are more consistent with a continuum of gene maturation than with a sharp division between "real genes" and prediction noise.

### Relationship to existing work

Our findings provide an independent computational line of support for the de novo gene birth model proposed by the Claverie group (Legendre et al. 2018; Jeudy et al. 2019), arrived at through a fundamentally different methodology — single-genome compositional analysis rather than multi-strain comparative genomics. The convergence of these independent approaches substantially strengthens the de novo gene birth hypothesis.

Our poly-A/T compositional analysis and Legendre et al.'s RNA secondary structure analysis may detect different aspects of the same termination system, operating at the sequence composition level and the RNA folding level respectively. The hairpin null model results (2.0 pp excess over null) suggest that simple palindrome detection captures only a fraction of the termination signal that energy-based RNA folding methods detect, and that the poly-A/T and hairpin mechanisms are complementary rather than competing.

The homology landscape reveals three distinct ecological narratives within the non-self hit taxonomy. The NCLDV connection (24 hits, 25.3%) represents the ancestral viral gene pool shared within the Nucleocytoviricota. The marine ecological network (28 algae + 13 Amoebozoa = 41 hits, 43.2%) reflects gene exchange within a shared coastal aquatic habitat — notably involving two phylogenetically distant algal lineages (green and brown algae), strengthening the ecological rather than phylogenetic interpretation. The endosymbiont signal (11 insect + 11 bacteria = 22 hits, 23.2%) is dominated by ankyrin repeats and likely reflects a shared bacterial HGT reservoir: aphid genomes contain many ankyrins acquired from Wolbachia endosymbionts, making the Pandoravirus–aphid connection indirect rather than direct.

The product-level gradient from ankyrin (42.3% non-self hit rate) through F-box (8.7%) to Fascin and Atrophin (0%) provides a temporal dimension: domains with high hit rates were likely acquired from external sources and retain recognizable homologs, while domains with zero external hits have either diverged beyond recognition or originated de novo within the Pandoravirus lineage.

### Additional compositional forces

The intergenic space from which new genes emerge is shaped by multiple forces beyond the AT regulatory system described here. Poirot et al. (2019) documented the near-complete absence of the AGCT tetranucleotide in Pandoravirus genomes, which we independently confirmed as a 1,500-fold depletion relative to Markov-1 null expectations (z = −82). Whether this reflects a restriction-modification system, a novel DNA editing process, or extreme compositional pressure, it demonstrates that Pandoravirus intergenic sequence is subject to selective constraints beyond neutral evolution. Sun et al. (2015) showed that MITE transposable elements colonize *P. salinus* intergenic space, raising questions about possible interactions between mobile elements and the distributed regulatory signals we describe. Whether MITEs disrupt, co-opt, or are recruited alongside AT-rich regulatory elements remains to be determined.

### Limitations

This study has several important limitations. All findings derive from a single genome (*P. salinus*); generalizability to other Pandoravirus strains and to other Nucleocytoviricota awaits confirmation. The transcriptomic validation used GFF annotations informed by RNA-seq data rather than raw read-level TSS/TTS coordinates, which are not publicly available; direct comparison of predicted regulatory signals against experimentally determined transcript boundaries would provide stronger validation. ESM-2 embeddings are a structural proxy, not actual 3D structure predictions; validation with AlphaFold or experimental structure determination would strengthen the structural maturation findings. The hypothetical gene subsample (n = 299) represents approximately one-third of all hypothetical genes; complete coverage would further refine the orphan rate estimates. Compositional methods reveal evolutionary history but cannot determine gene expression, functionality, or essentiality — the 92.5% ORFan fraction remains functionally uncharacterized. The CAI reference used general *A. castellanii* codon frequencies rather than a highly-expressed gene subset; CAI values are valid for relative comparisons but not absolute measures of translational optimization. Finally, the analysis pipeline was developed with AI tool assistance (see AI disclosure); while all results were independently verified, the analytical design reflects a human-AI collaborative process that may have introduced unrecognized biases in the choice of methods or interpretation of results.

### Predictions

Our model generates four falsifiable predictions. First, comparative genomic analysis of *P. dulcis* and *P. neocaledonia* should reveal proto-genes at different maturation stages, with some genes present in one strain corresponding to intergenic regions in others (partially confirmed by Jeudy et al. 2019 for *P. celtis*/*P. quercus*). Second, the strand-aware positional signal is consistent with promoter/terminator function; direct validation against experimentally determined TSS/TTS sites awaits reanalysis of primary transcriptomic data. Third, the 23 ORFan singletons that cluster structurally with known domain families (Fascin, ankyrin, kinase, Ubiquitin) should show structural similarity when analyzed with AlphaFold or experimentally determined. Fourth, we predict that poly-A/T runs in intergenic regions adjacent to hypothetical genes will show promoter or terminator activity in a reporter assay, even when the adjacent ORF is not conserved across Pandoravirus strains — providing a direct experimental test of the distributed regulatory infrastructure model.

---

## Use of AI Tools

Artificial intelligence tools (Anthropic's Claude) were used to assist with code development, data-processing workflow design, analytical planning, and editorial refinement of manuscript text. All analyses, statistical methods, interpretations, and final manuscript content were reviewed and validated by the author.

## Author Contributions

E.B. conceived the study, designed and implemented the analytical pipeline, performed all computational analyses, and wrote the manuscript.

## Data Availability

All analysis scripts, the PostgreSQL database schema, intermediate results, and raw BLAST output are available at [GitHub URL] and archived with a permanent DOI at [Zenodo DOI]. The *P. salinus* genome sequence is available under GenBank accession NC_022098.1. RNA-seq-informed GFF annotations were obtained from Legendre et al. (2018) via the IGS GBrowse server (http://www.igs.cnrs-mrs.fr/pandoraviruses/). A statistics audit file mapping every manuscript-cited number to its source computation is included in the repository.

---

## References

Aherfi S, Andreani J, Baptiste E, Oumessoum A, Dornas FP, Andrade ACDSP, Chafik A, Abrahao J, Colson P, La Scola B, Raoult D. 2018. A large open pangenome and a small core genome for giant pandoraviruses. *Frontiers in Microbiology* 9:1486.

Bisio H, Legendre M, Giry C, Philippe N, Alempic JM, Jeudy S, Abergel C. 2023. Evolution of giant pandoravirus revealed by CRISPR/Cas9. *Nature Communications* 14:428.

Jeudy S, Bertaux L, Alempic JM, Abergel C, Claverie JM. 2019. *Pandoravirus celtis* illustrates the microevolution processes at work in the giant *Pandoraviridae* genomes. *Frontiers in Microbiology* 10:430.

Legendre M, Fabre E, Poirot O, Jeudy S, Lartigue A, Alempic JM, Beucher L, Philippe N, Bertaux L, Christo-Foroux E, Labadie K, Byrne D, Léopold S, Abergel C, Claverie JM. 2018. Diversity and evolution of the emerging Pandoraviridae family. *Nature Communications* 9:2285.

Lin Z, Akin H, Rao R, Hie B, Zhu Z, Lu W, Smetanin N, Verkuil R, Kabeli O, Shmueli Y, dos Santos Costa A, Fazel-Zarandi M, Sercu T, Candido S, Rives A. 2023. Evolutionary-scale prediction of atomic-level protein structure with a language model. *Science* 379:1123–1130.

Philippe N, Legendre M, Doutre G, Couté Y, Poirot O, Lescot M, Arslan D, Seltzer V, Bertaux L, Bruber C, Criscuolo A, Jeudy S, Labadie K, Sterpenich A, Abergel C, Claverie JM. 2013. *Pandoraviruses*: amoeba viruses with genomes up to 2.5 Mb reaching that of parasitic eukaryotes. *Science* 341:281–286.

Poirot O, Jeudy S, Abergel C, Claverie JM. 2019. A Puzzling Anomaly in the 4-Mer Composition of the Giant Pandoravirus Genomes Reveals a Stringent New Evolutionary Selection Process. *Journal of Virology* 93(23):e01206-19.

Sun C, Feschotte C, Wu Z, Mueller RL. 2015. DNA transposons have colonized the genome of the giant virus *Pandoravirus salinus*. *BMC Biology* 13:38.

---

## Figure Legends

**Figure 1. Five-tier proto-gene continuum and maturation ordering.** (A) GC content across five confidence tiers, showing monotonic decline from T1 (0.650) to T5 (0.547), with the intergenic average (0.539) indicated by a dashed line. (B) Gene length distribution across tiers. (C) Codon usage displacement of T4 and T5 proto-genes relative to random intergenic ORFs and genome-optimized genes, showing 10–15% displacement toward genome-preferred codons. (D) Schematic of the four-stage maturation model: regulatory signals (pre-existing) → protein structure → codon optimization → GC composition.

**Figure 2. Distributed AT regulatory system and transcriptomic validation.** (A) Meta-gene AT content profile showing +10.3 pp spike at gene starts and +9.8 pp spike at gene ends. (B) Strand-aware classification of poly-A/T runs: 74% in correct regulatory position (binomial p < 10⁻⁶). (C) Comparison of AT boundary signal strength across Lambda, Mimivirus, and *P. salinus*. (D) Two-population test: boundary-proximal and dispersed runs show identical composition (AT ratio p = 0.105). (E) Proto-gene regulatory signals: T4 and T5 upstream AT content indistinguishable from annotated genes. (F) Hairpin null model: observed (11.6%) vs null (9.7%), 2.0 pp excess.

**Figure 3. BLASTp homology landscape.** (A) Taxonomic distribution of non-self top BLAST hits (n = 95), showing eukaryote-dominated pattern with marine algae (29.5%) and giant viruses (25.3%) as the largest categories. (B) Percent identity distribution of non-self hits (median 42.3%). (C) Product-level non-self hit rates, showing gradient from ankyrin (42.3%) through F-box (8.7%) to Pandoravirus-specific domains (0%). (D) Spatial distribution of homology across the genome, with non-self hits colored by taxonomy and native zone (1.3–2.0 Mb) highlighted. (E) Identity distribution by taxonomic category.

**Figure 4. Maturation gradient: annotated vs hypothetical homology.** (A) Non-self hit rates for annotated (18.2%, n = 521) versus hypothetical (3.0%, n = 299) genes, with Wilson score confidence intervals. Fisher's exact test OR = 7.19 (95% CI: 3.5–13.6, p = 7.79 × 10⁻¹²). The single highest-identity hypothetical hit (*Acanthamoeba castellanii*, 82.4% identity) is annotated as a likely recent host-to-virus HGT event. (B) Three-level maturation gradient: intergenic (0%) → hypothetical (3.0%) → annotated (18.2%). (C) E-value sensitivity analysis showing orphan rate stability across thresholds (81.8% at 10⁻⁵ to 90.4% at 10⁻⁵⁰).

---

## Supplementary Materials

**Table S1.** E-value sensitivity analysis for annotated gene non-self hit rates.

**Table S2.** Complete taxonomy of all 95 non-self BLAST hits with organism, identity, coverage, and genomic position.

**Table S3.** All 9 hypothetical gene non-self hits with taxonomy, identity, and native zone membership.

**Figure S1.** Full 7-panel Module 11 BLASTp analysis figure (panels A–G).

**Figure S2.** Module 10 corrected analysis 6-panel figure.
