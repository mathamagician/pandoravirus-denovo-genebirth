# Computational evidence for a distributed AT-based regulatory mechanism enabling de novo gene birth in *Pandoravirus salinus*

**Eddie Bradford**

Bradford Genomics, San Diego, CA, USA

---

## Abstract

*Pandoravirus salinus* possesses the largest known viral genome (2.47 Mb) with an unprecedented 92.5% of genes classified as ORFans, proteins with no detectable homologs outside the Pandoraviridae (Philippe et al. 2013). While de novo gene creation from non-coding sequence has been proposed as the primary diversifying mechanism of this lineage, the genomic infrastructure enabling this process remains poorly understood. Here we present a 12-module computational framework combining compositional analysis, structural prediction, transcriptomic validation, and genome-wide homology analysis to characterize the mechanism underlying de novo gene birth in *P. salinus*, with comparative validation across four additional Nucleocytoviricota spanning four families and GC contents from 28% to 64%.

We identify a five-tier proto-gene continuum exhibiting a monotonic gradient from established genes (mean GC 0.650, 1,452 bp) through proto-genes approaching intergenic composition (GC 0.547, 216 bp), with codon optimization detectable even in the most immature candidates (p = 0.0005). Critically, we discover that poly-A/T regulatory signals, functioning as promoters and terminators, pre-exist throughout intergenic space, with proto-genes exhibiting mature regulatory signatures indistinguishable from established genes (p = 0.14). This distributed AT-based regulatory system produces the strongest gene boundary signals observed across all tested genomes (+10.3 and +9.8 percentage points AT enrichment at gene starts and ends, respectively), despite operating against the genome's 62% GC composition.

Strand-aware analysis of poly-A/T positioning relative to experimentally informed gene annotations confirms that 74% of boundary-associated runs occupy the correct regulatory position (binomial p < 10⁻⁶). BLASTp analysis reveals that 81.8% of annotated genes and 97.0% of hypothetical genes are Pandoravirus-specific orphans, with annotated genes 7.2-fold more likely to retain external homologs than hypothetical genes (Fisher's exact test, OR = 7.19, 95% CI: 3.5–13.6, p = 7.79 × 10⁻¹²), a gradient that is robust across E-value thresholds (81.8%–90.4% from 10⁻⁵ to 10⁻⁵⁰).

Using an adaptive GC-normalized gene family detection method that corrects for base-composition bias in k-mer frequency space, we find that 64.1% of *P. salinus* genes are compositional singletons, consistent with de novo gene creation as the dominant growth mechanism. Comparative analysis across five NCLDVs demonstrates that several *P. salinus* patterns are broadly conserved: near-perfect strand balance in four of five genomes (48.8–51.3% plus-strand), translational selection toward host codon usage, and singleton-dominated proteomes. The direction of codon maturation, confirmed by formal Mann-Whitney U tests (GC3 and CAI shifts significant in all five genomes, all p < 0.05), reverses in AT-rich Mimivirus (28% GC), consistent with the de novo birth model across compositionally diverse genomes. *Pandoravirus dulcis* provides independent same-genus replication of the key *P. salinus* findings, and the AGCT tetranucleotide ranks as the single most depleted of all 256 4-mers (rank 1/256) in both Pandoravirus species and in *Mollivirus*, while enriched in *Mimivirus* (rank 241/256) and *Pithovirus* (rank 185/256), implicating a lineage-specific sequence constraint in the Pandoraviridae.

These findings establish a four-stage maturation model, regulatory signals (pre-existing) → protein structure → codon optimization → GC composition, providing independent computational confirmation of the de novo gene birth hypothesis via a fundamentally different methodology than the Claverie group's comparative genomic approach.

**Keywords:** de novo gene birth, ORFan genes, Pandoraviridae, regulatory signals, proto-gene continuum, genome evolution, Nucleocytoviricota

---

## Significance Statement

Giant viruses with genomes rivaling those of bacteria create most of their genes from scratch rather than acquiring them from other organisms, but how non-coding DNA gains the ability to function as a gene has remained unclear. We show that the largest known virus, *Pandoravirus salinus*, maintains a system of regulatory signals distributed throughout its non-coding DNA, providing the machinery for gene expression before genes even exist. This resolves a key puzzle, new genes can arise anywhere in the genome because the infrastructure to activate them is already everywhere, and reveals a specific evolutionary sequence from non-coding DNA to functional gene: first regulatory signals, then protein structure, then codon optimization, and finally bulk nucleotide composition. Comparative analysis across four additional giant viruses confirms that singleton-dominated proteomes and translational selection are broadly conserved features of the Nucleocytoviricota, while the extreme AGCT depletion and compositional uniformity are most pronounced in the Pandoraviridae.

---

## Introduction

The discovery of *Pandoravirus salinus* revealed the largest viral genome known to science, 2.47 megabases encoding approximately 1,430 predicted proteins, of which 92.5% have no detectable homologs outside the Pandoraviridae (Philippe et al. 2013). Our independent BLASTp analysis yields a consistent estimate of 91.4% (95% CI: 89.6%–93.0%; extrapolated from a subsample of 820 queried genes; see Results), confirming this extreme ORFan fraction by an orthogonal method. This unprecedented ORFan fraction poses a fundamental evolutionary question: where do these genes come from?

The Claverie group has provided compelling evidence that de novo gene creation from non-coding sequence is the primary diversifying mechanism of the Pandoraviridae. Legendre et al. (2018) conducted strand-specific RNA-seq across four Pandoravirus strains, demonstrating that 83–87% of genome sequence is transcribed, that palindromic hairpin structures mark approximately 70% of transcript 3' ends, and that the Pandoravirus pangenome is open, each new strain adds roughly 60 previously unseen proteins. Jeudy et al. (2019) provided direct evidence by comparing the nearly identical genomes of *P. celtis* and *P. quercus* (96.7% identity), identifying genes in one strain that correspond to intergenic regions in the other, genes captured in the act of being born. Aherfi et al. (2018) confirmed the open pangenome with a small conserved core, and Bisio et al. (2023) demonstrated through CRISPR/Cas9 knockout that essential genes concentrate at the 5' end while ORFans accumulate at the 3' end of *P. neocaledonia*.

These studies establish *what* is happening, genes are emerging from non-coding sequence, but leave a mechanistic gap: *how* does intergenic space acquire the regulatory and coding capacity to produce functional genes? If new genes require promoters and terminators to be expressed, how do they arise in intergenic regions distant from existing regulatory elements? And if gene birth is a gradual process, what is the sequence of events from non-coding DNA to functional protein?

Additional compositional analyses have revealed that Pandoravirus genomes are subject to unusual sequence constraints. Poirot et al. (2019) documented the near-complete absence of the AGCT tetranucleotide in some strains, pointing to a novel DNA editing or selection process. Sun et al. (2015) identified miniature inverted-repeat transposable elements (MITEs) colonizing the intergenic space of *P. salinus*, demonstrating that multiple forces shape the non-coding landscape from which new genes emerge.

Here we present a multi-method computational framework that addresses the mechanistic gap. Through a 12-module analytical pipeline encompassing genome architecture, null model calibration, codon usage analysis, gene prediction validation, evolutionary classification, gene family detection, protein structure prediction, regulatory signal discovery, transcriptomic validation, and genome-wide homology analysis, we provide four novel contributions. First, we establish that *P. salinus* is compositionally highly uniform despite being the largest genome tested, with ORFan genes indistinguishable from annotated genes by every codon metric applied, arguing against recent patchwork horizontal gene transfer. Second, we identify a quantitative five-tier proto-gene continuum with monotonic gradients across all measured features, supported by an adaptive-threshold gene family analysis and protein structure prediction, providing a maturation metric for gene birth. Third, we discover a distributed AT-based regulatory system that pre-loads intergenic regions with promoter and terminator signals, resolving the mechanistic gap of how genes arise far from existing regulatory elements. Fourth, we establish an asymmetric maturation ordering, regulatory signals pre-exist, followed by structural features, codon optimization, and finally GC composition, that defines the temporal sequence of de novo gene birth.

Comparative analysis of four additional NCLDV genomes (*Pandoravirus dulcis*, *Acanthamoeba polyphaga* mimivirus, *Pithovirus sibericum*, and *Mollivirus sibericum*) provides cross-lineage context, confirming which *P. salinus* findings are broadly conserved across the Nucleocytoviricota and which are lineage-specific to the Pandoraviridae.

Our approach differs fundamentally from previous work: all mechanistic findings derive from single-genome computational analysis of *P. salinus* without requiring multi-strain comparison, providing an independent line of evidence that complements the comparative genomics approach of the Claverie group.

---

## Results

### Genome-wide compositional landscape

To establish the baseline against which gene birth signatures are evaluated, we characterized the compositional architecture of *P. salinus* using a four-genome comparison ladder spanning three orders of magnitude in genome size (Table 1). *P. salinus* emerged as a compositional outlier within this set: it has the lowest coding density (67.5% vs 88–100% for the other three genomes), the largest intergenic fraction (32.5%), and by far the strongest GC gap between coding (64.6%) and intergenic (53.9%) regions at 10.7 percentage points, compared to 7.0 pp for Lambda and 4.1 pp for Mimivirus. This GC gap proved to be the most analytically productive architectural feature, as it reflects the AT-enriched regulatory signals concentrated in intergenic space.

**Table 1. Genome architecture across the four-genome comparison ladder.**

| Metric | PhiX174 | Lambda | Mimivirus | *P. salinus* |
|--------|---------|--------|-----------|-------------|
| Genome size (bp) | 5,386 | 48,502 | 1,181,549 | 2,473,870 |
| GC content | 44.7% | 49.9% | 28.0% | 61.7% |
| Coding fraction | 100.0% | 88.2% | 90.4% | 67.5% |
| Coding GC | 44.7% | 51.1% | 28.7% | 64.6% |
| Intergenic GC | N/A | 44.1% | 24.6% | 53.9% |
| GC gap (coding − IG) | N/A | 7.0 pp | 4.1 pp | **10.7 pp** |
| Avg gene length (bp) | 490 | 665 | 1,058 | 1,168 |
| Strand ratio (+/total) | 1.000 | 0.644 | 0.488 | 0.504 |

Null model calibration using Markov-1 models preserving dinucleotide frequencies revealed two sequence features of particular note. First, poly-A/T hexamers (AAAAAA/TTTTTT) are 15-fold enriched relative to null expectations (z > +100), confirming that homopolymer runs are genuine sequence features maintained by selection, not base composition artifacts. Second, the AGCT tetranucleotide is depleted 1,500-fold relative to Markov-1 expectations (z = −82), consistent with the extreme compositional constraint reported by Poirot et al. (2019). The effect-size filtering methodology proved essential: in a genome of 2.47 Mb, 92.6% of 4-mers pass z > 3 by statistical power alone, but adding an effect-size threshold (observed/expected < 0.5 or > 2.0) reduces this to 10.5%, demonstrating that 88.6% of nominally significant k-mer signals lack meaningful biological effect sizes.

Codon usage analysis revealed that *P. salinus* exhibits strong translational selection (neutrality plot R² = 0.045, N = 1,430 genes; the low R² is itself biologically informative, it indicates that GC at constrained positions, GC12, is largely invariant with changes in wobble GC, GC3, consistent with strong purifying selection at non-wobble sites; regression slope β = 0.118, 95% CI: 0.089–0.146, p = 5.9 × 10⁻¹⁶; the slope << 1.0 expected under strict neutrality confirms translational selection) with codons better adapted to the *Acanthamoeba* host translation machinery than Mimivirus (mean per-gene CAI 0.787 ± 0.032 for *P. salinus* vs 0.657 ± 0.027 for Mimivirus; Mann-Whitney U = 1,392,180, p effectively zero by double-precision overflow; mean difference 0.129, 95% bootstrap CI: 0.127–0.132). ORFan genes were compositionally indistinguishable from annotated genes across all codon metrics: GC3 difference +0.015, ENC difference −1.03, CAI difference +0.008, and upstream 6-mer cosine similarity 0.9824 (Supplementary Table S1). The compositional homogeneity indicates either ancient residency with full equilibration or de novo origin within the genome, not recent patchwork horizontal gene transfer.

Ab initio gene prediction with Prodigal (Hyatt et al. 2010) identified 1,685 CDS in *P. salinus*, confirming 93.4% of GenBank annotations (≥ 80% reciprocal overlap) while predicting 347 novel genes not in GenBank (mean length 380 bp, mean GC 57.0%). Note that Prodigal's hexamer log-likelihood scoring captures a different signal than our RSCU distance measure: Prodigal encodes translation initiation preference via hexamer composition at the start site, while RSCU distance measures genome-wide synonymous codon bias, the two statistics are complementary rather than redundant. The native zone (1.3–2.0 Mb) showed a higher confirmation rate (96.0% vs 92.4%) and lower hypothetical fraction (55.2% vs 66.2%) than the rest of the genome, it is better annotated, not a gene desert.

GMM clustering on five per-gene features (GC, GC3, GC12, ENC, CAI) identified three compositional clusters: a main native cluster (53.8%, GC 0.651), a high-GC cluster (29.4%, GC 0.654), and a variant cluster (16.8%, GC 0.613). These clusters primarily reflect quantitative GC variation within the genome rather than distinct evolutionary origins, and the 'NATIVE,' 'HIGH-GC,' and 'VARIANT' designations should be considered provisional pending multi-genome validation. Fifty-one extreme low-GC outliers (GC < 0.567), all hypothetical, were concentrated in the variant cluster. Despite this local variation, *P. salinus* exhibited a 5 kb sliding-window GC coefficient of variation of 3.9%, among the lowest of all tested genomes despite being the largest (Lambda: 14.2%, Mimivirus: 8.2%; Supplementary Table S8). This compositional uniformity is one of the strongest lines of evidence against recent patchwork assembly from diverse sources.

Gene family analysis using an adaptive GC-normalized threshold (see Materials and Methods) identified 12 composition-based gene groups in *P. salinus*, encompassing 36% of genes; the remaining 64.1% are compositional singletons with no detectable similar partners. Of the 12 groups, the most biologically interpretable are small families of repeat-domain proteins (Fascin-like: 12 members; ankyrin: 8, 4, and 2 members) with no tandem clustering, members are dispersed across the genome. The largest group (478 genes, predominantly hypothetical proteins of similar GC content and length) shows genuine sequence-level coherence above null expectation: subgroup analysis using composite 6-mer Jaccard and sliding-window percent identity reveals large effect sizes above singleton controls for both the 208 domain-containing members (Cohen's d = 2.62, p = 4.9 × 10⁻⁹¹) and the 270 hypothetical protein members (d = 2.55, p = 7.0 × 10⁻⁸⁸). These effect sizes substantially exceed the near-zero d = 0.06 observed for random gene pairs under the composition-matched null, indicating that the group reflects genuine distant paralogy rather than compositional convergence. The subgroup also contains a preponderance of repeat-domain proteins (ankyrin, F-box, MORN, DHFR-like: 208 of 478), consistent with a MITE-derived expansion of short repetitive elements as documented by Sun et al. (2015). Alignment-based methods (MMseqs2, CD-HIT, or HMM profiles) are needed to delineate specific sub-lineages within this group; for now it is best treated as a compositionally coherent set of likely distant paralogs. Null calibration confirmed that the diffuse pairwise similarity cloud of the full proteome (mean 0.109 vs Mimivirus 0.069) is fully explained by *P. salinus*'s skewed amino acid composition driven by 62% GC (null mean 0.105; Cohen's d = 0.06), while the repeat-domain groups themselves rise above a well-characterized null background (Supplementary Table S3).

For cross-genome comparison, it is important to note that fixed compositional similarity thresholds perform inconsistently across genomes with different GC content. At 28% GC, the null distribution mean for random 3-mer pairs exceeds 0.95, meaning the commonly used 0.95 threshold would classify nearly all Mimivirus gene pairs as similar, a clear artifact. The adaptive threshold corrects this bias and enables valid comparison. Under adaptive thresholds, singleton rates across the five NCLDV genomes tested range from 64% to 91% (Supplementary Table S4), indicating that de novo gene creation, not duplication, is the primary growth mechanism across this lineage sample.

### A five-tier proto-gene continuum

Classification of all 2,292 predicted ORFs (1,430 GenBank + 862 additional Prodigal predictions) into five confidence tiers revealed a monotonic gradient across every measured compositional feature (Table 2). GC content declined continuously from 0.650 in mature annotated genes (T1) to 0.547 in the most immature proto-gene candidates (T5), approaching the intergenic average of 0.539. Gene length followed the same trajectory (1,452 bp → 216 bp), as did GC3 (0.750 → 0.672 for the tiers where wobble position could be computed) and CAI (0.792 → 0.762).

**Table 2. Five-tier proto-gene continuum.** Spearman ρ = −1.0 for both GC content and gene length across five tiers (N = 5, exact one-tailed p = 0.0083 for each; Fisher combined p = 0.00073 for the joint test). A perfect rank correlation with N = 5 is statistically uncommon by chance (1 in 120 permutations), but the biological significance lies in the monotonicity across independent measurements (GC, length, GC3, CAI) rather than in any single p-value.

| Tier | Description | n | Mean length (bp) | Mean GC | GC3 | CAI |
|------|------------|---|------------------|---------|-----|-----|
| T1 | Annotated, Prodigal-confirmed | 516 | 1,452 | 0.650 | 0.750 | 0.792 |
| T2 | Hypothetical, Prodigal-confirmed | 784 | 1,082 | 0.648 | 0.745 | 0.787 |
| T3 | GenBank-only (Prodigal missed) | 130 | 555 | 0.615 | 0.672 | 0.762 |
| T4 | Prodigal high-confidence only | 335 | 501 | 0.560 |, |, |
| T5 | Prodigal low-confidence only | 527 | 216 | 0.547 |, |, |

The gradient did not merely reflect annotation quality. Codon usage analysis comparing T4 and T5 predictions against length-matched random intergenic ORFs revealed statistically significant genome-directed optimization even in the least mature candidates: T4 codon usage was significantly closer to the genome-preferred pattern than random (p = 0.020), with T5 showing an even stronger signal (p = 0.0005). The effect size was modest but precisely measured: T4 predictions showed a mean normalized RSCU distance of 0.74 toward the genome reference versus 0.76 for length-matched random ORFs (6.7% displacement toward the genome-optimized baseline of 0.46); T5 showed 1.01 versus 1.03 (3.5% displacement), consistent with early-stage but real genome-directed optimization in recently emerged ORFs. Prodigal's gene-finding algorithm uses hexamer log-likelihood scoring rather than explicit codon preference optimization (Hyatt et al. 2010), so the RSCU-based signal we detect represents a partially independent measure of coding character that must originate from the sequences themselves rather than from the prediction algorithm's scoring criteria.

Proto-genes showed no spatial clustering near established genes (T4 p = 0.99; T5 p = 1.0), ruling out regulatory co-option or operon extension as mechanisms. Instead, proto-gene candidates were distributed uniformly across the genome, consistent with a model in which gene birth can occur anywhere in the intergenic landscape.

ESM-2 protein language model embeddings for 170 representative proteins spanning all compositional clusters and proto-gene tiers revealed 9 structural families. Null calibration against composition-matched random proteins confirmed that real structural families are significantly tighter than expected by chance (within-family similarity 0.533 vs null 0.368; Cohen's d = 1.45), though 41% of PC1 variance correlates with GC content, indicating a partial compositional confound. Twenty-three ORFan singletons clustered structurally with annotated domain-containing proteins, suggesting deeper homology relationships invisible to sequence-based methods (Supplementary Table S2). Eight of nine structural families were cross-tier, indicating that structural relationships cut across compositional categories.

*P. dulcis* replicates the *P. salinus* five-tier pattern at the genus level. Comparative analysis reveals essentially identical compositional architecture: GC CV of 3.6% vs 3.9%, similar coding density (67.9% vs 67.5%), and positive GC3 shift of annotated versus hypothetical genes (+0.036 vs +0.016), all consistent with the same de novo maturation pathway operating independently in a second Pandoravirus species. This replication transforms the proto-gene continuum from a single-species observation to a plausible Pandoraviridae-level property.

### A distributed AT-based regulatory system

The most striking finding of the compositional analysis was the discovery of an AT-based regulatory system operating against the genome's high GC content. After filtering to intergenic runs only (see Transcriptomic validation below), *P. salinus* exhibited gene boundary AT signals dramatically stronger than any other genome in our comparison set: +10.3 percentage points AT enrichment at gene starts and +9.8 percentage points at gene ends, compared to approximately +1 percentage point for Lambda and +1–2 percentage points for Mimivirus. This signal runs counter to the genome's 62% GC composition, requiring active regulatory selection to maintain. The full cross-genome comparison of AT boundary signals is provided in Supplementary Table S10.

Poly-A and poly-T homopolymer runs (≥ 5 bp) numbered 6,955 and 7,010 respectively, with 82.8–83.9% located in intergenic regions, a pattern unique to *P. salinus* in this comparison set (Lambda: 12–24%; Mimivirus: 16–17%). These runs displayed strand-asymmetric positioning consistent with regulatory function: poly-A runs were enriched upstream of gene starts (41.5% within 200 bp), while poly-T runs were enriched downstream of gene ends (40.8%).

Divergent intergenic regions (flanked by genes transcribed outward in both directions) were 17% longer than convergent regions (582 bp vs 498 bp; p < 0.0001), consistent with a requirement for two promoters to be accommodated. However, divergent regions were less dense in poly-A runs (6.67 vs 7.69 runs/kb; p = 0.019), indicating that the extra length represents spacer DNA between two normal-density promoter regions.

The critical finding for the de novo gene birth model was that proto-gene candidates exhibited regulatory signals indistinguishable from established genes. T4 upstream AT content was 45.0% (p = 0.14 vs annotated 45.9%) and T5 was 46.7% (p = 0.25), neither significantly different from mature genes. This means that the regulatory infrastructure for transcription pre-exists in intergenic space before ORFs emerge, addressing the mechanistic gap of how spatially random proto-genes can be expressed: they are born into a landscape already equipped with promoters and terminators.

### Transcriptomic validation

To validate the computationally predicted regulatory signals, we compared poly-A/T run positions against experimentally informed gene annotations from Legendre et al. (2018). After filtering to 11,640 intergenic runs and applying strand-aware classification, we found that 74% of boundary-associated runs occupied the correct regulatory position, poly-A upstream on the sense strand or poly-T upstream on the antisense strand for promoters, and the converse for terminators (binomial p < 10⁻⁶).

This strand-aware analysis identified 2,332 runs in promoter-consistent positions and 1,227 in terminator-consistent positions. The asymmetry (1.9:1 promoter-to-terminator ratio) may reflect the greater complexity of transcription initiation relative to termination, or may indicate that a subset of termination events uses the hairpin mechanism described by Legendre et al. (2018) rather than poly-T signals.

Meta-gene AT content profiles, averaging AT content across aligned gene boundaries, revealed sharp spikes at transcript boundaries: +10.3 percentage points at gene starts and +9.8 percentage points at gene ends. These values exceed the initial unfiltered estimates (+7.0/+9.5 pp from Module 9) because the intergenic filter removed coding-region runs that diluted the signal.

To assess whether palindromic hairpin structures, proposed by Legendre et al. (2018) as the primary termination signal, could account for the poly-A/T signal, we implemented a hairpin null model with strict parameters (stem ≥ 6 bp, exact complementarity). The observed hairpin fraction (11.6%) exceeded the null model expectation (9.7%) by only 2.0 percentage points (z = 3.6, statistically significant but biologically small). The 70% hairpin frequency reported by Legendre et al. (2018) likely reflects energy-based RNA folding analysis (e.g., RNAfold/mfold) rather than simple palindrome detection. Importantly, the poly-A/T signal stands independently of hairpin structures.

The two-population test comparing boundary-proximal (n = 4,331) and dispersed intergenic (n = 7,309) poly-A/T runs found no meaningful compositional difference: mean AT ratio difference of 1.6 percentage points (p = 0.105) and no significant length difference (p = 0.093). A 1.6 pp effect on a background mean of 38.3% AT corresponds to a negligible effect size (Cohen's d < 0.05), confirming that boundary-proximal and dispersed intergenic runs are effectively compositionally identical. This supports the latent regulatory substrate model: dispersed runs throughout intergenic space are compositionally equivalent to those actively functioning at gene boundaries. The regulatory potential is distributed, not concentrated, any region of the genome could, in principle, support transcription of a newly emerged ORF.

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
| Other virus | 1 | 1.1 |, |

The marine algae category represents an ecological grouping, not a phylogenetic one: green algae (Chlorophyta/Prasinophyceae; *Micromonas*, *Bathycoccus*, *Ostreococcus*) and brown algae (Phaeophyceae; *Ectocarpus*, *Scytosiphon*) are as distantly related as animals and fungi, yet both share marine coastal habitat with *P. salinus*. The 11 insect hits mapped exclusively to ankyrin and F-box domains in aphid genomes (*Acyrthosiphon pisum*, *Aphis* spp.); these likely reflect convergent domain expansion mediated by Wolbachia endosymbionts rather than direct Pandoravirus-to-insect horizontal transfer.

Non-self hits were spatially enriched in the native zone (1.3–2.0 Mb): 45 of 95 hits fell within this region (OR = 1.94, Fisher's exact p = 0.004), consistent with the native zone containing evolutionarily older genes that retain broader phylogenetic connections. Product-level analysis revealed a gradient from clearly acquired domains (ankyrin repeat: 42.3% non-self hit rate) through diverged domains (F-box: 8.7%; MORN repeat: 4.5%) to Pandoravirus-specific domains with zero external homologs (Fascin, Atrophin, DHFR, Ring finger, BTB: 0%), providing direct evidence for different evolutionary ages of domain acquisition or divergence (Supplementary Table S3).

Of 299 successfully queried hypothetical proteins, only 9 (3.0%) had non-self hits, compared to 95 of 521 (18.2%) for annotated proteins, a 7.2-fold enrichment (Fisher's exact test, OR = 7.19, 95% CI: 3.5–13.6, p = 7.79 × 10⁻¹²; Table 4). This result was robust across E-value thresholds: the annotated orphan rate remained above 80% even at the most stringent threshold tested (81.8% at 10⁻⁵, 83.1% at 10⁻¹⁰, 84.1% at 10⁻²⁰, 90.4% at 10⁻⁵⁰; Supplementary Table S7).

**Table 4. Annotated vs hypothetical homology comparison.**

| Gene category | n | Non-self hit rate | Interpretation |
|--------------|---|------------------|----------------|
| Annotated (recognized domains) | 521 | 18.2% | Oldest genes; some retain detectable homology |
| Hypothetical (ORFs, no domains) | 299 | 3.0% | Younger genes; almost entirely orphans |
| Intergenic regions | ~1,900 IGRs | 0% (by definition) | Not yet genes; the raw substrate |

The 9 hypothetical non-self hits were also spatially enriched in the native zone: 7 of 9 (77.8%) fell within 1.3–2.0 Mb versus 67 of 290 (23.1%) without hits (OR = 11.7, Haldane-corrected OR = 9.98, 95% CI: 2.32–42.83, p = 0.0010). The 9 hypothetical non-self hits displayed a taxonomic profile distinct from annotated hits: 6 giant viruses (66.7%), 2 Amoebozoa (22.2%), and 1 bacterium (11.1%), with zero marine algae and zero insects. Hypothetical genes with detectable external homology connect exclusively to phylogenetically proximate sources (the NCLDV gene pool and the amoebozoan host lineage), while more distant connections (marine algae, insects, diverse bacteria) appear only in annotated genes. This phylogenetic distance gradient, younger genes retain similarity only to nearest relatives, is an additional prediction of the maturation model confirmed by the data.

### Cross-NCLDV comparative analysis

To determine which *P. salinus* findings reflect broadly conserved NCLDV properties versus Pandoraviridae-specific features, we applied the compositional pipeline (Modules 2–6 and 7-v2) to four additional giant virus genomes spanning GC contents from 28% to 64% and four distinct families (Supplementary Table S6).

Three patterns are consistent across all five genomes tested: (1) near-perfect strand balance in four of five genomes (48.8–51.3% plus-strand for Mimivirus, *P. salinus*, *P. dulcis*, and Mollivirus), with *Pithovirus sibericum* showing a modest minus-strand bias at 44.3%, approximately 2.5 standard deviations from 50%, possibly reflecting replication-associated strand asymmetry or a distinct gene orientation constraint; (2) translational selection toward host codon usage, with neutrality slopes ranging from 0.067 to 0.173 for four of the five genomes (slopes well below the neutral expectation of 1.0); and (3) singleton-dominated proteomes, with adaptive-threshold singleton rates of 64–91% across the full GC range, indicating that de novo gene creation, not duplication, is the primary growth mechanism across this lineage sample. These broadly conserved features are consistent with all five genomes independently deploying the same fundamental strategy for gene origination.

The direction of codon maturation reverses with genome base composition, providing a key cross-validation of the maturation model, and this reversal is confirmed by formal Mann-Whitney U tests. In all four GC-rich genomes (*P. salinus*, *P. dulcis*, *Mollivirus*, *Pithovirus*), annotated genes have higher GC3 and higher CAI than hypothetical genes, consistent with maturation driving codons toward the high-GC genome-preferred usage. In Mimivirus (28% GC), annotated genes have *lower* GC3 (Δ = −0.024, Mann-Whitney p = 1.09 × 10⁻⁹, Cohen's d = −0.439) and *lower* CAI (Δ = −0.012, p = 2.30 × 10⁻¹³, d = −0.467) than hypothetical genes, consistent with maturation toward AT-preferred codons in a low-GC genome. This reversal is precisely what the de novo gene birth model predicts: regardless of GC content, older (more mature) genes are better adapted to their host's translational machinery. Among the GC-rich genomes, the strongest maturation effects are observed in *Mollivirus* (GC3: p = 1.01 × 10⁻⁹, d = 0.693; CAI: p = 2.01 × 10⁻¹⁰, d = 0.703) and *P. dulcis* (GC3: p = 8.3 × 10⁻⁸, d = 0.469; CAI: p = 7.5 × 10⁻¹⁰, d = 0.467). *P. salinus* shows the smallest effect (GC3: p = 0.027, d = 0.227; CAI: p = 7.0 × 10⁻⁵, d = 0.246), and *Pithovirus* shows significant GC3 shift (p = 4.7 × 10⁻⁴, d = 0.354) but no significant CAI shift (p = 0.850, d = 0.027), consistent with its higher compositional heterogeneity (GC CV 13.5%) and more complex evolutionary history. The reversal would not be expected under a model of recent horizontal acquisition, which would produce foreign composition signatures regardless of maturation state. Full statistics for all five genomes are provided in Supplementary Table S6.

*Pandoravirus dulcis* provides the strongest available validation for the *P. salinus* findings. As an independent species from the same genus, it replicates the key compositional patterns: GC CV of 3.6% vs *P. salinus* 3.9%, positive GC3 maturation direction (+0.036 vs +0.016), similar coding density (67.9% vs 67.5%), and AGCT near/complete absence (0 copies in 1.9 Mb vs 4 copies in *P. salinus* 2.5 Mb). Regulatory signal analysis of *P. dulcis* confirms that the distributed AT system operates at the genus level: 86% of poly-A/T runs are intergenic (matching the 82.8–83.9% range in *P. salinus*), gene boundary AT enrichment is +7.7 pp at start positions, and 61.8% of boundary-associated runs occupy the correct strand-aware regulatory position, all fully consistent with *P. salinus* values. The compositional architecture is so similar that several claims in this paper, compositional homogeneity, positive GC3 maturation direction, singleton-dominated proteome, distributed AT regulatory system, are now supported at the genus level rather than the single-species level.

Exhaustive BLASTp analysis of all *P. dulcis* proteins (319 annotated valid queries; 741 hypothetical valid queries) confirms genus-level homology patterns consistent with *P. salinus*: annotated orphan rate 94.7% and hypothetical orphan rate 98.4%, with annotated genes significantly more likely to retain non-self homologs than hypothetical genes (OR = 3.42, Haldane-corrected OR = 3.38, p = 1.51 × 10⁻³). Notably, the *P. dulcis* annotated orphan rate (94.7%) is significantly higher than the *P. salinus* rate (81.8%; Fisher's exact test, OR = 0.253, Haldane-corrected OR = 0.259, p = 2.86 × 10⁻⁸), indicating that *P. dulcis* annotated genes are more lineage-specific. This difference is unlikely to reflect methodological artefact, both analyses were exhaustive and used identical protocols, and is consistent with *P. dulcis* having a proportionally larger complement of recently emerged, lineage-restricted genes, or with its annotated genes having diverged further from their ancestral homologs. The hypothetical orphan rates are indistinguishable between species (98.4% vs 97.0%), confirming that the upper bound on orphan fraction is a conserved genus-level property.

*Pithovirus sibericum* (36% GC) stands out as compositionally the most heterogeneous genome in the set (GC CV 13.5%, compared to 3–8% for the other four), consistent with a more mosaic evolutionary history potentially involving larger HGT contributions. This compositional heterogeneity makes Pithovirus an attractive contrast system for future HGT-versus-de-novo-birth analysis.

One anomaly requires explicit note. *Mollivirus sibericum* shows a negative neutrality slope (−0.146), unique among the five genomes. Unlike the other four genomes where low positive slopes indicate that constrained codon positions are nearly invariant (translational selection), a negative slope means GC3 and GC12 move in opposite directions, GC3 decreases as GC12 increases. This could reflect amino acid-level selection opposing nucleotide composition pressure, an unusual compositional dynamic in a 60% GC genome with 80.7% hypothetical proteins, or an artifact of poor annotation in a recently described species with limited functional characterization. *Mollivirus* should be treated as an anomalous case in neutrality analyses rather than grouped with the four genomes showing positive translational selection. Notably, *Mollivirus* does show significant GC3 and CAI maturation shifts (d = 0.693 and 0.703, respectively), indicating that the annotated-vs-hypothetical maturation signal is robust even in a genome whose overall neutrality plot is anomalous.

---

## Discussion

### A mechanistic model for de novo gene birth

The results presented here converge on a coherent model for de novo gene birth in *P. salinus* that can be articulated as a three-step logical chain. First, poly-A/T runs with regulatory positioning exist throughout the intergenic space, as demonstrated by the 82.8–83.9% intergenic fraction, the strand-asymmetric boundary enrichment, and the +10.3/+9.8 pp meta-gene AT spikes confirmed by strand-aware validation. Second, boundary-proximal and dispersed intergenic runs are compositionally identical (1.6 pp AT ratio difference, negligible effect size, p = 0.105), indicating that the regulatory substrate is distributed uniformly rather than concentrated at existing gene boundaries. Third, proto-genes emerge into these pre-existing regulatory contexts with mature transcriptional signals already in place (T4 p = 0.14, T5 p = 0.25), acquire coding features through a gradual maturation process measurable as the five-tier continuum, and progressively accumulate external homologs as they age (3.0% → 18.2% non-self hit rate).

This model predicts a specific maturation ordering. Regulatory signals pre-exist as a feature of intergenic sequence composition, representing the earliest (indeed, prerequisite) stage. Protein structural features emerge next, proto-genes were not structurally separable from established domain families in embedding space (nearest-neighbor similarity test, p = 0.90), consistent with early acquisition of protein-level organization even among the most immature candidates. Codon optimization follows, detectable as a 3.5–6.7% displacement from random toward genome-preferred usage (T4 RSCU distance 0.74 vs random 0.76, p = 0.020; T5 distance 1.01 vs random 1.03, p = 0.0005). GC composition shifts last, as the slowest and most passive process, with proto-genes (GC 0.547) only gradually approaching the coding average (GC 0.650) from the intergenic baseline (GC 0.539).

It is important to note that this ordering is inferred from a cross-sectional snapshot of genes at different maturation stages, not from longitudinal tracking of individual genes through time. Whether this sequence is a necessary ordering (i.e., regulatory signals must precede structural features, which must precede codon optimization) or merely the typical ordering observed in this genome cannot be determined from these data alone. A gene could in principle acquire partial codon optimization before fully establishing structural features, if the two processes are not strictly coupled. The cross-sectional consistency of the gradient across all tiers is consistent with the proposed ordering, but confirming necessity versus typicality requires comparative genomic tracking across multiple Pandoravirus strains at different evolutionary distances.

### Evidence against alternative genome growth models

Multiple independent lines of evidence converge against the two standard alternative models for genome expansion. Against recent patchwork horizontal gene transfer: *P. salinus* exhibits very low 5 kb sliding-window GC coefficient of variation (3.9%) despite being the largest genome tested, replicating in *P. dulcis* (GC CV 3.6%), and confirming this as a Pandoraviridae-level compositional property; ORFan genes are compositionally indistinguishable from annotated genes across all codon metrics (GC3 difference +0.015, CAI difference +0.008); upstream regulatory signatures are shared (6-mer cosine similarity 0.9824); GMM clustering finds no dramatically foreign compositional outliers; and the BLASTp taxonomy shows an eukaryote-dominated pattern consistent with shared ecology rather than multi-kingdom HGT.

Against recent gene duplication: adaptive-threshold gene family analysis identifies 64.1% of *P. salinus* genes as compositional singletons with no similar partners; the repeat-domain families involve domain-level expansions (ankyrin, Fascin-like), not whole-gene duplications, and family members are dispersed across the genome with no tandem clustering; and the diffuse pairwise similarity cloud (mean 0.109) is fully explained by a null model matched for base composition (null mean 0.105; Cohen's d = 0.06). The largest detected group (478 genes) shows genuine sequence-level coherence above background (domain subgroup d = 2.62, hypothetical subgroup d = 2.55), consistent with distant paralogy, likely a MITE-mediated repeat expansion (Sun et al. 2015), but this accounts for only ~33% of genes; the dominant mode of genome expansion remains non-duplicative at 64.1% singletons. The same singleton-dominated architecture is observed across all five NCLDVs tested (64–91% with adaptive thresholds), consistent with a shared evolutionary strategy in which genome expansion proceeds primarily through de novo gene origination rather than duplication.

An alternative interpretation is that the T4 and T5 proto-gene tiers reflect annotation artifacts or false-positive Prodigal predictions rather than genuine evolutionary intermediates. Several observations argue against this. First, both T4 and T5 exhibit codon usage significantly shifted toward genome-preferred patterns relative to length-matched random intergenic ORFs, with the strongest signal retained even in T5 (p = 0.0005), indicating sequence-level structure inconsistent with random noise. Second, representative proto-genes are not structurally isolated but cluster near established domain families in embedding space, suggesting continuity in protein-level organization rather than arbitrary ORF detection. Third, proto-genes show no spatial clustering near established genes, arguing against simple operon extension or boundary spillover. We acknowledge a partial circularity in the tier classification: annotation practice assigns lower confidence to shorter, lower-GC sequences, and Prodigal preferentially detects shorter ORFs in AT-richer regions, the same properties that define the lower tiers. A definitive separation of annotation bias from biological maturation signal requires comparative genomic validation across multiple Pandoravirus strains.

### Relationship to existing work

Our findings provide independent computational confirmation of the de novo gene birth model proposed by Legendre et al. (2018) and Jeudy et al. (2019), arrived at through a fundamentally different methodology, single-genome compositional analysis rather than multi-strain comparative genomics. The convergence of these independent approaches substantially strengthens the de novo gene birth hypothesis.

Our poly-A/T compositional analysis and the RNA secondary structure analysis of Legendre et al. (2018) may detect different aspects of the same termination system. The hairpin null model results (2.0 pp excess over null) suggest that simple palindrome detection captures only a fraction of the termination signal that energy-based RNA folding methods detect, and that the poly-A/T and hairpin mechanisms are complementary rather than competing.

The homology landscape reveals three distinct ecological narratives within the non-self hit taxonomy. The NCLDV connection (24 hits, 25.3%) represents the ancestral viral gene pool shared within the Nucleocytoviricota. The marine ecological network (28 algae + 13 Amoebozoa = 41 hits, 43.2%) reflects gene exchange within a shared coastal aquatic habitat; the prevalence of marine algal sequences among NCLDV homologs is consistent with broad ecological exchange between giant viruses and their algal hosts documented by Moniruzzaman et al. (2020) and with the taxonomic framework described by Aylward et al. (2021). The endosymbiont signal (11 insect + 11 bacteria = 22 hits, 23.2%) is dominated by ankyrin repeats and likely reflects a shared bacterial HGT reservoir mediated by Wolbachia endosymbionts. The product-level gradient from ankyrin (42.3% non-self hit rate) through F-box (8.7%) to Fascin and Atrophin (0%) provides a temporal dimension for domain acquisition or divergence.

### Cross-NCLDV comparative validation

The five-genome comparative analysis reveals a consistent picture across the sampled NCLDV diversity, with important caveats about the sample size and family representation. Three features appear broadly conserved: strand balance near 50% (with the exception of *Pithovirus*, which shows a modest minus-strand bias at 44.3%), translational selection toward host codon usage, and singleton-dominated proteomes. The consistency of these features across four NCLDV families and a 36-percentage-point range in GC content suggests they reflect fundamental constraints of the Nucleocytoviricota infection strategy, though confirmation across a wider genome sample is needed before universal claims are warranted.

The GC3 maturation reversal between AT-rich Mimivirus and the four GC-rich genomes is conceptually the most important comparative finding. It reframes the maturation model in a genome-composition-independent way: maturation drives codons toward *host-adapted* usage, not toward any fixed GC content. Formal Mann-Whitney U tests confirm the reversal is statistically robust, with large effect sizes in both directions (Mimivirus d = −0.44 for GC3, −0.47 for CAI; *Mollivirus* d = +0.69 and +0.70; *P. dulcis* d = +0.47 and +0.47). The old Mimivirus singleton rate of 9.5%, produced by a fixed cosine threshold applied below the random noise floor for AT-rich genomes, was not merely inaccurate, it was meaningfully wrong in the other direction. The corrected rate of 88.9% makes the cross-genome picture coherent: all five genomes show singleton-dominated proteomes with adaptive thresholds. This is the key finding that the method correction enables.

*P. dulcis* is the single most valuable addition to the analysis from the perspective of the *P. salinus* paper. It converts the central claim from "this *P. salinus* pattern is interesting" to "*P. salinus* and *P. dulcis* independently show the same pattern", a qualitatively different level of evidence. For the Pandoravirus-specific claims (compositional homogeneity, AGCT extreme depletion, positive GC3 maturation direction, distributed AT regulatory system), *P. dulcis* provides the nearest-possible validation: a distinct species from the same genus, independently analyzed, showing the same results.

### Additional compositional forces

The intergenic space from which new genes emerge is shaped by multiple forces beyond the AT regulatory system described here. The AGCT depletion we confirmed in *P. salinus* (1,500-fold relative to Markov-1 expectations, z = −82) proves to be part of a broader gradient across giant viruses, but with a Pandoraviridae-specific extreme. Ranking AGCT among all 256 4-mers by observed/expected ratio reveals a clear gradient: AGCT is modestly enriched in Mimivirus (rank 241/256; O/E = 1.27) and Pithovirus (rank 185/256; O/E = 1.11); the single most depleted 4-mer of all 256 in Mollivirus (rank 1/256; O/E = 0.11); and essentially or completely absent in both Pandoravirus species, rank 1/256 in *P. salinus* (O/E = 0.00066, 4 copies across 2.5 Mb) and rank 1/256 in *P. dulcis* (O/E = 0.00, 0 copies across 1.9 Mb). The rank-1/256 status in three genomes of varying GC content confirms that AGCT depletion is not merely a GC-driven compositional consequence: Mollivirus (60% GC) and the Pandoraviridae (62–64% GC) both show extreme AGCT depletion, while Pithovirus (36% GC) and Mimivirus (28% GC) show modest enrichment. The comparison with Mollivirus at similar GC content (60%) but 259 AGCT copies versus *P. salinus*'s 4 copies demonstrates that the Pandoraviridae depletion substantially exceeds what GC content alone predicts, implicating a lineage-specific mechanism.

Whether this pattern reflects a restriction-modification system, a novel DNA editing process, or extreme compositional pressure, it demonstrates that Pandoravirus intergenic sequence is subject to selective constraints beyond neutral evolution. The restriction-modification hypothesis remains speculative at this stage; the data establish the pattern and its lineage specificity, not the mechanism. The MITE colonization documented by Sun et al. (2015) raises questions about possible interactions between mobile elements and the distributed regulatory signals we describe. Whether MITEs seeded the AT regulatory landscape during colonization of *P. salinus* intergenic space, or whether pre-existing AT-rich regulatory sequences facilitated MITE insertion, or whether both processes co-occurred remains to be resolved.

### Limitations

This study has several important limitations. The mechanistic findings derive primarily from *P. salinus*, and while the five-genome comparison provides useful context, five genomes from four families remain a small sample of the Nucleocytoviricota. Patterns described as "broadly conserved across sampled NCLDV diversity" should not be interpreted as universal without confirmation across a wider range of families. *Mollivirus* in particular, with an anomalous negative neutrality slope and 80.7% hypothetical proteins reflecting limited functional annotation, should be treated with additional caution in cross-genome comparisons.

The gene family analysis warrants several explicit caveats. The 3-mer cosine similarity method, even with the adaptive GC-normalized threshold, is a compositional clustering tool, not a sequence homology detector. Clusters reflect compositional similarity, shared k-mer frequencies, which may correlate with but does not confirm true paralogy. The 478-gene group, now supported by subgroup sequence similarity analysis (d = 2.55–2.62), most likely reflects distant paralogy, possibly MITE-related repeat expansion, rather than compositional noise alone; however, alignment-based validation (MMseqs2, CD-HIT, HMM profiles) is required to delineate specific paralogy sub-lineages within this group and to determine how many distinct ancestral sequences gave rise to its members. Benchmarking the adaptive 3-mer method against established family detection algorithms on at least Mimivirus and *P. salinus* remains a planned but incomplete step. Sensitivity analysis at the 95th and 99.5th percentile thresholds confirms that qualitative conclusions, high singleton rates, no massive whole-gene duplication, are robust to threshold choice (Supplementary Table S5). Cross-genome comparisons of singleton rates should also account for the fact that annotation quality varies substantially: Mollivirus at 80.7% hypothetical proteins and Pithovirus at 72.8% are less well-characterized than Mimivirus at 54.0%.

To clarify the scope of the comparative analysis: Modules 2–7 (compositional architecture, null models, codon usage, gene prediction, GMM, gene families) were applied to all five NCLDV genomes. Module 8 (ESM-2 protein structure) was applied to *P. salinus* only for the primary mechanistic findings. Module 9 (regulatory signal analysis) was applied to *P. salinus* as the primary genome of study, with comparative analysis additionally performed on *P. dulcis* to validate genus-level generalization. Module 11 (BLASTp homology analysis) was completed exhaustively for both *P. salinus* and *P. dulcis*. Comparative claims therefore reflect different levels of evidence: all five-genome claims are cross-genome Module 2–7 results, while protein structure findings are *P. salinus*-specific and regulatory signal findings are *P. salinus* primary with *P. dulcis* validation.

The transcriptomic validation used GFF annotations informed by RNA-seq data rather than raw read-level TSS/TTS coordinates, which are not publicly available; direct comparison of predicted regulatory signals against experimentally determined transcript boundaries would provide stronger validation. ESM-2 embeddings are a structural proxy, not actual 3D structure predictions. The hypothetical gene subsample (n = 299) represents approximately one-third of all hypothetical genes; the 91.4% overall orphan estimate is extrapolated from this subsample and should be interpreted accordingly. Compositional methods reveal evolutionary history but cannot determine gene expression, functionality, or essentiality. The CAI reference used general *A. castellanii* codon frequencies rather than a highly-expressed gene subset; all five viruses analyzed here infect *Acanthamoeba*, making this a consistent reference across the dataset, but CAI values are valid for relative comparisons and not absolute measures of translational optimization. Multiple statistical tests were performed across 12 analytical modules; we distinguish confirmatory tests (pre-specified as direct tests of central hypotheses: proto-gene codon usage comparison, strand-aware 74% binomial test, annotated-vs-hypothetical Fisher's exact test, and two-population AT run test) from exploratory comparisons (k-mer enrichment profiles, cross-genome metrics, spatial subgroup analyses). All within-genome gene comparisons assume independence of individual gene statistics, but genes share evolutionary history; p-values for within-genome comparisons should be interpreted as measures of consistency and effect magnitude rather than calibrated false-positive rates.

Finally, the analysis pipeline was developed with AI tool assistance (see Use of AI Tools); while all results were independently verified, the analytical design reflects a human-AI collaborative process.

### Predictions

Our model generates four falsifiable predictions. First, comparative genomic analysis across Pandoravirus strains should reveal proto-genes at different maturation stages, with some genes in one strain corresponding to intergenic regions in others (partially confirmed by Jeudy et al. 2019 for *P. celtis*/*P. quercus*; our *P. dulcis* comparative analysis confirms the same compositional maturation signatures, consistent with this prediction). Second, the strand-aware positional signal is consistent with promoter/terminator function; direct validation against experimentally determined TSS/TTS sites awaits reanalysis of primary transcriptomic data. Third, the 23 ORFan singletons that cluster structurally with known domain families (Fascin, ankyrin, kinase, Ubiquitin) should show structural similarity when analyzed with AlphaFold or experimentally determined. Fourth, we predict that poly-A/T runs in intergenic regions adjacent to hypothetical genes will show promoter or terminator activity in a reporter assay, even when the adjacent ORF is not conserved across Pandoravirus strains, providing a direct experimental test of the distributed regulatory infrastructure model.

---

## Materials and Methods

### Genome and database

The complete genome of *Pandoravirus salinus* (GenBank accession NC_022098.1; 2,473,870 bp, 1,430 annotated CDS) was loaded into a PostgreSQL 15 database alongside three comparison genomes spanning the dsDNA virus size spectrum: bacteriophage PhiX174 (NC_001422.1; 5,386 bp), bacteriophage Lambda (NC_001416.1; 48,502 bp), and *Acanthamoeba polyphaga* mimivirus (NC_014649.1; 1,181,549 bp).

For cross-NCLDV comparative validation, three additional giant virus genomes were analyzed using the same pipeline: *Pandoravirus dulcis* (NC_021858.1; 1,908,524 bp, 1,070 annotated CDS), *Pithovirus sibericum* (NC_023423.1; 610,033 bp, 467 annotated CDS), and *Mollivirus sibericum* (NC_027867.1; 651,523 bp, 523 annotated CDS). Together, the five NCLDV genomes span GC contents from 28.0% to 63.7%, genome sizes from 610 kb to 2.5 Mb, and four distinct families (Mimiviridae, Pandoraviridae, Pithoviridae, Molliviridae). The database schema comprises 22 tables including position-indexed nucleotide sequences, gene annotations, intergenic regions, codon usage statistics, and derived analytical tables. All analyses used Python 3.11 with BioPython 1.86, NumPy, SciPy, scikit-learn, and matplotlib.

### Compositional and structural analysis pipeline (Modules 2–8)

A 12-module analytical pipeline was applied to characterize the compositional and structural landscape of *P. salinus*, with Modules 2–6 and 7-v2 additionally applied to all five NCLDV genomes for comparative analysis. Module 8 (protein structure) was applied to *P. salinus* as the primary genome of study. Module 9 (regulatory signal analysis) was applied to *P. salinus* with comparative analysis additionally performed for *P. dulcis*.

**Genome architecture (Module 2).** Eight metrics were computed genome-wide and per-gene: GC content, coding density, strand bias, intergenic region length, Shannon entropy, gene length distribution, and 5 kb sliding-window GC profiles. The GC coefficient of variation reported throughout this paper (3.9% for *P. salinus*, 3.6% for *P. dulcis*) refers specifically to this 5 kb sliding-window metric, which captures spatial GC uniformity across the genome rather than per-gene variation (per-gene GC CV is higher at ~6.1%). All metrics were compared across the four-genome ladder and, for the compositional metrics, across all five NCLDVs.

**Null model calibration (Module 3).** Markov-1 null models preserving dinucleotide frequencies were used to assess k-mer enrichment and depletion (k = 4, 6, 8). Effect-size filtering (observed/expected < 0.5 or > 2.0) was applied alongside z-score significance to reduce false positives in large genomes, where 88.6% of z > 3 signals lacked meaningful effect sizes. To confirm that mononucleotide-based null models are adequate for M7 threshold calibration, we additionally implemented dinucleotide-preserving (Markov-1) random sequence generation and compared the resulting null thresholds against mononucleotide-based thresholds; the two approaches differed by Δ = +0.0001 in the 99th-percentile threshold, confirming that mononucleotide null calibration is adequate for gene family detection.

**Codon usage analysis (Module 4).** Per-gene relative synonymous codon usage (RSCU), GC3 (wobble position), GC12 (constrained positions), effective number of codons (ENC), and codon adaptation index (CAI) were computed. CAI was calculated using *Acanthamoeba castellanii* codon frequencies from the Kazusa database; all five NCLDV genomes analyzed here infect *Acanthamoeba*, making this a consistent host reference across the comparative analysis. Annotated versus hypothetical gene comparisons used Mann-Whitney U tests with Cohen's d effect sizes to test for compositional signatures indicative of foreign acquisition or maturation stage.

**Gene prediction validation (Module 5).** Ab initio gene prediction with Prodigal (Hyatt et al. 2010) was compared against GenBank annotations using a reciprocal overlap criterion: a Prodigal prediction was considered to confirm a GenBank annotation if (overlap length / GenBank annotation length) ≥ 0.80 AND (overlap length / Prodigal prediction length) ≥ 0.80; both conditions must be met simultaneously. Novel Prodigal-only predictions were characterized by length, GC content, and genomic position. Note that Prodigal's hexamer log-likelihood scoring (capturing translation initiation site preference) measures a different signal than our per-gene RSCU distance (measuring genome-wide synonymous codon bias); the two statistics are complementary rather than redundant.

**Evolutionary origin classification (Module 6).** Gaussian Mixture Model (GMM) clustering on five per-gene features (GC, GC3, GC12, ENC, CAI) identified compositional clusters. These clusters should be interpreted as quantitative GC variation within the genome rather than confirmed distinct evolutionary origins; the biological interpretation of GMM classes as evolutionary categories requires orthogonal validation. The native zone (1.3–2.0 Mb), defined independently by this clustering, was used in subsequent spatial analyses.

**Gene family detection, adaptive threshold (Module 7-v2).** All-vs-all protein 3-mer frequency cosine similarity was computed for each genome. Fixed cosine similarity thresholds are susceptible to base-composition bias: in AT-rich genomes, the k-mer frequency space is compressed, inflating pairwise similarity among unrelated sequences. At 28% GC (Mimivirus), the null distribution mean for random 3-mer pairs exceeds 0.95, meaning a 0.95 fixed threshold classifies essentially all gene pairs as similar, a catastrophic false-positive artifact. To correct this, we developed an adaptive GC-normalized threshold: for each genome, 5,000 pairs of random sequences are generated matching the genome's base composition and gene length distribution, pairwise 3-mer cosine similarities are computed on these random pairs, and the 99th percentile of the resulting null distribution is used as the clustering threshold. This threshold adapts automatically to any genome's GC content, ensuring a consistent ~1% false-positive rate under the null. Sensitivity analysis confirms that qualitative conclusions are robust at the 95th and 99.5th percentile thresholds (Supplementary Table S5). Note that this method produces composition-based groups, not sequence homology families; clusters reflect shared k-mer composition, which may co-occur with but does not confirm true paralogy. For the 478-gene *P. salinus* group, we additionally performed subgroup sequence similarity analysis using composite 6-mer Jaccard and sliding-window percent identity metrics, comparing within-group pairs to singleton controls (300 sampled pairs each) to assess whether compositional grouping reflects genuine sequence-level coherence. Validation against alignment-based family detection methods (MMseqs2, CD-HIT) is planned. Hierarchical clustering via connected components identified groups of ≥2 genes above the adaptive threshold; genes with no above-threshold partners are classified as compositional singletons.

**Proto-gene tier classification (Modules 7a–7c).** All 2,292 predicted ORFs were classified into five confidence tiers (T1–T5) based on prediction agreement between GenBank annotations and Prodigal. For the disambiguation codon usage test, the RSCU Euclidean distance to the genome-preferred codon vector (mean RSCU across all T1+T2 genes) was computed for each T4 and T5 candidate and compared against 10,000 length-matched random intergenic ORFs using one-sided Mann-Whitney U tests.

**Protein structure prediction (Module 8).** ESM-2 (650M parameters; Lin et al. 2023) embeddings were generated for 170 representative *P. salinus* proteins spanning all tiers and compositional clusters. PCA reduction (top 20 components, 85.8% variance) followed by hierarchical clustering identified structural families. Null calibration confirmed structural signal genuineness (within-family similarity 0.533 vs null 0.368; Cohen's d = 1.45), though 41% of PC1 variance correlates with GC content, indicating a partial compositional confound.

### Regulatory signal analysis (Module 9)

Seven complementary analyses characterized the regulatory landscape of *P. salinus* and *P. dulcis*: (i) genome-wide poly-A/T run inventory (runs ≥ 5 bp), (ii) meta-gene AT content profiles across gene boundaries, (iii) upstream k-mer enrichment analysis, (iv) intergenic grammar comparison between convergent and divergent regions, (v) strand-asymmetric boundary positioning, (vi) ORFan versus annotated regulatory signature comparison, and (vii) proto-gene regulatory context assessment. The same pipeline was applied to *P. dulcis* for comparative genus-level analysis.

### Transcriptomic validation (Module 10)

Computationally predicted regulatory signals were compared against experimentally informed gene annotations from Legendre et al. (2018). GFF annotations incorporating RNA-seq evidence were used in lieu of raw reads, which are not deposited in public repositories.

Three methodological corrections were applied relative to an initial analysis: (i) an intergenic filter removed 2,325 poly-A/T runs (16.6%) overlapping coding sequences, retaining 11,640 intergenic runs; (ii) a strand-aware boundary test classified each run as promoter-positioned or terminator-positioned, replacing a strand-naive proximity test; (iii) a hairpin null model using strict parameters (stem ≥ 6 bp, exact complementarity) quantified the contribution of palindromic sequences to termination signals.

The two-population test compared AT composition and length distributions between boundary-proximal runs (within 200 bp of a gene boundary; n = 4,331) and dispersed intergenic runs (n = 7,309) to test whether boundary-active and dispersed runs represent distinct populations or a uniform regulatory substrate. Both the AT ratio difference (1.6 pp on a background mean of 38.3%) and the length difference were non-significant and negligible in effect size, supporting the uniform substrate model.

### BLASTp homology analysis (Modules 11–12)

All 528 annotated (non-hypothetical) *P. salinus* protein sequences were searched against the NCBI non-redundant protein database (nr) using remote BLASTp via BioPython (E-value threshold 1e-5, top 10 hits per query). Of 528 queries, 521 returned valid results; 7 failed queries due to server timeouts were excluded. A subsample of 300 randomly selected hypothetical proteins (from 902 total; random seed 42 for initial 100, seed 43 for additional 200) was searched under identical parameters, yielding 299 successful queries.

Non-self hits were defined as hits to organisms other than Pandoravirus strains. Taxonomic classification used a curated organism-level lookup table covering all 50 unique organisms encountered in the dataset, supplemented by an ordered keyword fallback and NCBI Entrez API lookup for unrecognized organisms. During quality auditing, the initial automated classification was found to have misidentified 12 prasinophyte green algae as bacteria due to a substring match on the suffix '-monas'; this was corrected by implementing the curated lookup table approach.

The annotated versus hypothetical comparison used Fisher's exact test with Haldane-corrected Woolf method for the 95% confidence interval on the odds ratio, and Wilson score intervals for proportions. E-value sensitivity analysis re-filtered existing hits at thresholds of 10⁻¹⁰, 10⁻²⁰, and 10⁻⁵⁰ to assess robustness.

---

## Use of AI Tools

Artificial intelligence tools (Anthropic's Claude) were used to assist with code development, data-processing workflow design, analytical planning, and editorial refinement of manuscript text. All analyses, statistical methods, interpretations, and final manuscript content were reviewed and validated by the author.

---

## Author Contributions

E.B. conceived the study, designed and implemented the analytical pipeline, performed all computational analyses, and wrote the manuscript.

---

## Data Availability

All analysis scripts, the PostgreSQL database schema, intermediate results, and raw BLAST output are available at https://github.com/mathamagician/pandoravirus-denovo-genebirth and archived with a permanent DOI at https://doi.org/10.5281/zenodo.19347286. The *P. salinus* genome sequence is available under GenBank accession NC_022098.1. RNA-seq-informed GFF annotations were obtained from Legendre et al. (2018) via the IGS GBrowse server (http://www.igs.cnrs-mrs.fr/pandoraviruses/). A statistics audit file mapping every manuscript-cited number to its source computation is included in the repository.

---

## References

Aherfi S, Colson P, La Scola B, Raoult D. 2018. A large open pangenome and a small core genome for giant Pandoraviruses. *Front Microbiol* 9:1486.

Aylward FO, et al. 2021. Formal recognition of the Kingdom Bamfordvirae and two new orders to accommodate large and giant DNA viruses of the phylum Nucleocytoviricota. *PLOS Biol* 19:e3001430.

Bisio H, Legendre M, Giry C, Abergel C, Claverie JM. 2023. Functional characterization of a giant virus through gene deletion analysis. *Nat Commun* 14:428.

Hyatt D, Chen GL, Locascio PF, Land ML, Larimer FW, Hauser LJ. 2010. Prodigal: prokaryotic gene recognition and translation initiation site identification. *BMC Bioinformatics* 11:119.

Jeudy S, Bertaux L, Alempic JM, Abergel C, Claverie JM. 2019. *Pandoravirus celtis* illustrates the microevolution processes at work in the giant *Pandoraviridae* genomes. *Front Microbiol* 10:430.

Legendre M, et al. 2018. Diversity and evolution of the emerging Pandoraviridae family. *Nat Commun* 9:2285.

Lin Z, et al. 2023. Evolutionary-scale prediction of atomic-level protein structure with a language model. *Science* 379:1123–1130.

Moniruzzaman M, Weinheimer AR, Martinez-Gutierrez CA, Aylward FO. 2020. Widespread endogenization of giant viruses shapes genomes of green algae. *Nature* 588:141–145.

Philippe N, et al. 2013. *Pandoraviruses*: amoeba viruses with genomes up to 2.5 Mb reaching that of parasitic eukaryotes. *Science* 341:281–286.

Poirot O, Jeudy S, Abergel C, Claverie JM. 2019. A puzzling anomaly in the 4-mer composition of the giant Pandoravirus genomes reveals a stringent new evolutionary selection process. *J Virol* 93(23):e01206-19.

Sun C, Feschotte C, Wu Z, Mueller RL. 2015. DNA transposons have colonized the genome of the giant virus *Pandoravirus salinus*. *BMC Biol* 13:38.

---

## Figure Legends

**Figure 1. Five-tier proto-gene continuum and maturation ordering.** (A) GC content across five confidence tiers, showing monotonic decline from T1 (0.650) to T5 (0.547), with the intergenic average (0.539) indicated by a dashed line. (B) Gene length distribution across tiers (box plots showing median, interquartile range, and outliers). (C) Codon usage displacement of T4 and T5 proto-genes relative to random intergenic ORFs and genome-optimized annotated genes (T1+T2), showing 3.5–6.7% displacement toward genome-preferred codons. (D) Schematic of the four-stage maturation model: regulatory signals (pre-existing) → protein structure → codon optimization → GC composition.

Alt text: Four-panel figure showing the proto-gene continuum. Panel A is a bar chart of mean GC content declining from tier 1 (0.650) to tier 5 (0.547). Panel B shows box plots of gene length decreasing across tiers. Panel C is a bar chart showing codon usage displacement of proto-genes between random and genome-optimized levels. Panel D is a schematic diagram with four sequential arrows showing the maturation order.

**Figure 2. Distributed AT regulatory system and transcriptomic validation.** (A) Meta-gene AT content profile showing +10.3 pp spike at gene starts and +9.8 pp spike at gene ends. (B) Strand-aware classification of poly-A/T runs: 74% in correct regulatory position (binomial p < 10⁻⁶). (C) Comparison of AT boundary signal strength across Lambda, Mimivirus, and *P. salinus*. (D) Two-population test comparing boundary-proximal runs (within 200 bp of a gene boundary, n = 4,331) against dispersed intergenic runs (n = 7,309): mean AT ratio difference 1.6 pp (p = 0.105, Cohen's d < 0.05), supporting a uniform regulatory substrate. (E) Proto-gene regulatory signals: T4 and T5 upstream AT content indistinguishable from annotated genes. (F) Hairpin null model: observed (11.6%) vs null (9.7%), 2.0 pp excess.

Alt text: Six-panel figure showing the AT regulatory system. Panel A is a line plot of AT content across gene boundaries showing sharp upward spikes at start and end positions. Panel B is a stacked bar chart showing 74% of poly-A/T runs in correct strand-aware regulatory position. Panel C is a grouped bar chart comparing AT boundary signal strength across three genomes. Panel D shows overlapping density plots for boundary-proximal and dispersed runs with nearly identical distributions. Panel E is a bar chart comparing upstream AT content across gene tiers. Panel F is a bar chart comparing observed and null-model hairpin frequencies.

**Figure 3. BLASTp homology landscape.** (A) Taxonomic distribution of non-self top BLAST hits (n = 95), showing eukaryote-dominated pattern with marine algae (29.5%) and giant viruses (25.3%) as the largest categories. (B) Percent identity distribution of non-self hits (median 42.3%). (C) Product-level non-self hit rates, showing gradient from ankyrin (42.3%) through F-box (8.7%) to Pandoravirus-specific domains (Fascin, DHFR, Atrophin, Ring finger: 0%). (D) Spatial distribution of non-self homology hits across the genome, with points colored by taxonomy and the native zone (1.3–2.0 Mb, shaded) highlighted (45/95 hits, OR = 1.94, p = 0.004). (E) Identity distribution by taxonomic category.

Alt text: Five-panel figure showing the BLASTp homology landscape. Panel A is a pie chart of taxonomic categories of non-self hits. Panel B is a histogram of percent identity values. Panel C is a horizontal bar chart of product-level hit rates. Panel D is a scatter plot across the genome with a shaded native zone. Panel E shows box plots of identity distributions by taxonomic category.

**Figure 4. Maturation gradient: annotated vs hypothetical homology.** (A) Non-self hit rates for annotated (18.2%, n = 521) versus hypothetical (3.0%, n = 299) genes, with Wilson score confidence intervals. Fisher's exact test OR = 7.19 (95% CI: 3.5–13.6, p = 7.79 × 10⁻¹²). (B) Three-level maturation gradient: intergenic (0%) → hypothetical (3.0%) → annotated (18.2%). (C) E-value sensitivity analysis showing orphan rate stability across thresholds (81.8% at 10⁻⁵ to 90.4% at 10⁻⁵⁰).

Alt text: Three-panel figure showing the maturation gradient. Panel A is a bar chart comparing non-self hit rates with confidence intervals. Panel B is a step chart showing the three-level gradient. Panel C is a line plot showing orphan rates remaining above 80% across increasingly stringent E-value thresholds.

---

## Supplementary Materials

Supplementary Tables S1–S10 and Supplementary Figures S1–S9 are provided in a single supplementary PDF file.

**Supplementary Table S1.** Codon usage comparison between annotated and ORFan genes: RSCU, GC3, GC12, ENC, and CAI with statistical tests.

**Supplementary Table S2.** ESM-2 structural family assignments for 170 representative proteins, including tier, compositional cluster, and nearest-neighbor similarity.

**Supplementary Table S3.** Composition-based gene group detection results: all 12 groups with member counts, mean similarity, spatial distribution, biological annotation, and product-level non-self BLASTp hit rates.

**Supplementary Table S4.** Adaptive GC-normalized gene family detection across five NCLDV genomes: null distribution parameters, adaptive thresholds, group counts, and singleton rates.

| Genome | GC% | Null mean | Adaptive threshold (99th pctile) | Groups (≥2) | Singleton rate |
|--------|-----|-----------|----------------------------------|-------------|----------------|
| *A. polyphaga* mimivirus | 28.0% | 0.9593 | 0.9905 | 16 | 88.9% |
| *Pithovirus sibericum* | 35.8% | 0.9332 | 0.9827 | 4 | 90.4% |
| *Mollivirus sibericum* | 60.1% | 0.9296 | 0.9825 | 16 | 81.6% |
| *P. salinus* | 61.7% | 0.9467 | 0.9840 | 12 | 64.1% |
| *P. dulcis* | 63.7% | 0.9513 | 0.9863 | 5 | 68.7% |

**Supplementary Table S5.** Sensitivity analysis: singleton rates at 95th, 99th, and 99.5th percentile adaptive thresholds across all five NCLDV genomes.

**Supplementary Table S6.** Cross-NCLDV comparative metrics summary with formal Mann-Whitney U statistics for GC3/CAI maturation shifts.

| Metric | Mimivirus | *P. salinus* | *P. dulcis* | Pithovirus | Mollivirus |
|--------|-----------|------------|-----------|------------|------------|
| GC content | 28.0% | 61.7% | 63.7% | 35.8% | 60.1% |
| Genome length | 1.18 Mb | 2.47 Mb | 1.91 Mb | 610 kb | 651 kb |
| Genes | 979 | 1,430 | 1,070 | 467 | 523 |
| Orphan fraction | 54.0% | 63.1% | 69.8% | 72.8% | 80.7% |
| Coding density | 89.5% | 67.5% | 67.9% | 68.0% | 82.6% |
| Strand balance (+ strand) | 48.8% | 50.4% | 51.3% | 44.3%‡ | 49.1% |
| GC CV (5 kb windows, %) | 8.2 | 3.9 | 3.6 | 13.5 | 3.2 |
| GC3 shift (ann − hyp) | −0.024 | +0.016 | +0.036 | +0.022 | +0.052 |
| GC3 Mann-Whitney p | 1.09×10⁻⁹ | 0.027 | 8.3×10⁻⁸ | 4.7×10⁻⁴ | 1.01×10⁻⁹ |
| GC3 Cohen's d | −0.439 | 0.227 | 0.469 | 0.354 | 0.693 |
| CAI shift (ann − hyp) | −0.012 | +0.008 | +0.016 | +0.001 | +0.027 |
| CAI Mann-Whitney p | 2.30×10⁻¹³ | 7.0×10⁻⁵ | 7.5×10⁻¹⁰ | 0.850† | 2.01×10⁻¹⁰ |
| CAI Cohen's d | −0.467 | 0.246 | 0.467 | 0.027† | 0.703 |
| Neutrality slope | 0.067 | 0.114 | 0.142 | 0.173 | −0.146§ |
| AGCT O/E ratio | 1.27× | 0.0007× | 0× | 1.11× | 0.11× |
| AGCT rank (of 256 4-mers) | 241 | 1 | 1 | 185 | 1 |
| Singleton rate (adaptive) | 88.9% | 64.1% | 68.7% | 90.4% | 81.6% |

†Pithovirus CAI shift is not statistically significant (p = 0.850, d = 0.027); GC3 shift is significant (p = 4.7×10⁻⁴, d = 0.354).
‡Pithovirus shows a modest minus-strand bias (~2.5σ from 50%); all other four genomes are 48.8–51.3%.
§Mollivirus shows an anomalous negative neutrality slope; see Discussion.

**Supplementary Table S7.** E-value sensitivity analysis for annotated gene non-self hit rates.

**Supplementary Table S8.** Cross-genome compositional homogeneity: GC, GC3, GC12, ENC, and CAI coefficient of variation comparison.

**Supplementary Table S9.** Hypothesis testing matrix: evidentiary status of four genome growth models (HGT, gene duplication, mobile element expansion, de novo gene birth) across all 12 modules, with supporting/disfavoring classification and key statistics.

**Supplementary Table S10.** Gene boundary AT signals across genomes: AT content, start and end delta-AT values, and poly-A/T intergenic fraction for Lambda, Mimivirus, *P. salinus*, and *P. dulcis*.

**Supplementary Figure S1.** Genome architecture: sliding-window GC profiles, gene density maps, GC skew plots, and strand maps for all four genomes in the comparison ladder.

**Supplementary Figure S2.** GMM clustering of *P. salinus* genes: three compositional clusters with native zone overlay.

**Supplementary Figure S3.** Pairwise protein 3-mer cosine similarity distributions: *P. salinus* vs Mimivirus, with composition-matched null model overlay demonstrating that the diffuse similarity cloud (mean 0.109) is a compositional artifact (null mean 0.105; Cohen's d = 0.06).

**Supplementary Figure S4.** ESM-2 PCA embedding space showing structural family assignments, proto-gene positions, GC content overlay, and null calibration comparison (within-family similarity 0.533 vs null 0.368; Cohen's d = 1.45).

**Supplementary Figure S5.** Adaptive threshold null distributions for all five NCLDV genomes, showing the relationship between genome GC content and the null similarity mean. At 28% GC (Mimivirus), the null mean (0.9593) exceeds the fixed threshold of 0.95, demonstrating the failure mode of fixed-threshold methods in AT-rich genomes.

**Supplementary Figure S6.** Cross-NCLDV GC3 shift direction reversal: annotated vs hypothetical gene GC3 and CAI across five genomes, showing positive shift in four GC-rich genomes and negative shift in Mimivirus, with effect sizes (Cohen's d) and Mann-Whitney p-values.

**Supplementary Figure S7.** AGCT depletion gradient across five NCLDVs, with Markov-1 expected values and rank-among-256 annotation. The Pandoraviridae (*P. salinus*, *P. dulcis*) show depletion substantially beyond what GC content alone predicts, compared against *Mollivirus* at similar GC content.

**Supplementary Figure S8.** Full 7-panel Module 11 BLASTp analysis figure (panels A–G).

**Supplementary Figure S9.** Module 10 corrected analysis 6-panel figure.

---

