# Technical Walkthrough: Pandoravirus de novo Gene Birth Analysis
## For statistically sophisticated readers with limited genomics background

**Eddie Bradford — Pandoravirus salinus GBE Submission, Draft 6**

---

## Orientation: What is this paper actually about?

Before diving into methods, it helps to understand what we are trying to explain and why it is puzzling.

**The organism.** *Pandoravirus salinus* is a giant virus — a virus with a genome so large (2.47 megabases, roughly half the size of a minimal bacterial genome) that it encodes ~1,430 proteins. Most viruses encode a few dozen. This one is an outlier.

**The puzzle.** 92.5% of its proteins are "ORFans" — they have no detectable relatives anywhere in the known biological world. This is extraordinary. For comparison, even highly unusual bacteria typically have 80–90% of their genes with detectable homologs somewhere. Where did all these unique genes come from?

**The proposed answer.** The leading hypothesis is *de novo gene birth*: the virus creates new genes from scratch, out of what was previously non-coding "junk" DNA in between existing genes. A stretch of random-looking intergenic sequence gains the properties of a functional gene over evolutionary time.

**The mechanistic gap this paper addresses.** Even if you accept de novo gene birth as the answer, a puzzle remains: for a new gene to be useful, it needs to be *expressed* — the cell's machinery needs to be told "start reading here, stop reading there" (these signals are called promoters and terminators). How does a piece of junk DNA get expressed if it is far from any existing regulatory signal? This paper's central finding is that the answer is already everywhere: Pandoravirus maintains regulatory signals distributed throughout its entire non-coding genome, so any stretch of sequence that randomly develops coding potential is already sitting next to a promoter/terminator and can immediately be expressed.

---

## Biological vocabulary, translated for statisticians

| Biology term | Statistical/analytical equivalent |
|---|---|
| **Gene** | A coded sequence that produces a protein; the "unit of analysis" |
| **ORFan / orphan gene** | A gene with zero matches in external databases; a complete outlier in the homology space |
| **Intergenic region (IGR)** | The sequence *between* genes — non-coding DNA, the substrate from which new genes emerge |
| **Promoter / terminator** | Sequence signals that tell the cell where to start/stop reading a gene — like delimiters in a file |
| **Codon** | A 3-letter DNA code word that specifies one amino acid; there are 64 codons but only 20 amino acids, so most amino acids are encoded by multiple synonymous codons |
| **GC content** | The fraction of the DNA letters that are G or C (vs A or T) — the primary compositional axis |
| **Annotated gene** | A gene with a known function or recognized protein domain, validated by bioinformatics databases |
| **Hypothetical gene** | A gene-like region that has no recognized function — an ORF (open reading frame) with no known relatives |
| **Proto-gene** | A gene-like sequence at an intermediate stage of maturation — not fully gene-like, not fully random |
| **BLASTp** | A database search algorithm: for a given protein sequence, find the most similar sequences across all known biology. Returns hits ranked by E-value. |
| **NCLDV** | Nucleocytoviricota — the taxonomic group of giant viruses; Pandoravirus and Mimivirus are both NCLDVs |

---

## The data structure

Everything starts with the genome sequence of *P. salinus*: a single circular DNA molecule of 2,473,870 base pairs. This was stored in a PostgreSQL database with 22 tables tracking:

- Position-indexed sequence (every nucleotide, position-labeled)
- Gene annotations (start, end, strand, protein sequence, function label)
- Intergenic regions (the gaps between genes)
- Derived analytical tables (codon usage statistics, k-mer frequencies, etc.)

The genome was analyzed alongside three comparison genomes chosen to span the size spectrum: bacteriophage PhiX174 (5,386 bp), bacteriophage Lambda (48,502 bp), and Mimivirus (1,181,549 bp). For cross-lineage comparison, four additional giant virus genomes were added: *Pandoravirus dulcis*, *Pithovirus sibericum*, *Mollivirus sibericum*, and Mimivirus.

The analytical unit is usually the individual gene: 1,430 annotated CDS in *P. salinus*, giving a comfortable sample size for gene-level statistics. For sequence-level analyses, the unit is typically a k-mer (a sequence window of length k) counted across the whole genome.

---

## Module 1: Null model calibration — separating signal from base-composition noise

**Why this matters first.** *P. salinus* has 62% GC content — meaning 62% of its DNA letters are G or C. This extreme base composition creates a massive statistical nuisance: almost any sequence-level enrichment or depletion you might find could just be an artifact of the unusual letter frequencies, not a real biological signal. The very first analytical step was to build null models that bake in this base composition, so that downstream signals are genuine.

**Markov-1 null model.** We used a Markov-1 (first-order) model — a model that preserves observed *dinucleotide* frequencies (every pair of adjacent letters) rather than just individual letter frequencies. This is one step more realistic than a simple Bernoulli model of independent letters, because consecutive nucleotides in real genomes are not independent.

For each k-mer of interest, we computed: observed count vs expected count under the Markov-1 model, and a z-score. This immediately reveals two extreme patterns:

- **Poly-A/T hexamers (AAAAAA / TTTTTT): 15-fold enriched** (z > +100). These runs are dramatically over-represented relative to even a GC-matched null — they exist because selection maintains them, not because the base composition produces them by chance.

- **AGCT tetranucleotide: 1,500-fold depleted** (z = −82). This specific 4-letter sequence is almost completely absent from the genome. With only 4 copies across 2.47 million base pairs, this is one of the most extreme single-sequence depletion patterns in known biology.

**The effect-size filter.** Here is the key methodological contribution: in a genome of 2.47 Mb, statistical power is so high that 92.6% of all 4-mers pass z > 3 (α ≈ 0.001) by raw significance alone. But a z-score of 3.2 on a 3-fold-enriched 4-mer is not the same as a z-score of 82 on a 1,500-fold-depleted 4-mer. Filtering to |O/E| > 2 (at least a 2-fold deviation from null expectation) reduces the significant fraction from 92.6% to 10.5% — meaning 88.6% of nominally significant signals lack biologically meaningful effect sizes. The paper consistently applies this effect-size discipline throughout, distinguishing statistical significance from practical significance.

---

## Module 2: Codon usage analysis — measuring translational adaptation

**Background.** Every amino acid in a protein is specified by a 3-letter codon. Most amino acids have 2–6 synonymous codons that encode the exact same thing. Different organisms preferentially use different synonymous codons — this "codon bias" reflects adaptation to the organism's translation machinery (the machinery that reads DNA and makes proteins). A gene that is well-adapted to its host's translation machinery uses the host's preferred codons; a newly-born gene that hasn't been "optimized" by selection yet will use more random codon choices.

Three metrics were computed per gene:

**GC3 (wobble position GC content).** The third position of each codon is usually "silent" — changing it doesn't change the amino acid. This position is called the "wobble" position. The GC content at wobble positions (GC3) is a clean readout of compositional pressure on synonymous sites, relatively free from amino-acid-level selection.

**CAI (Codon Adaptation Index).** A 0-to-1 score measuring how well a gene's codon choices match the codon frequencies of the host organism (*Acanthamoeba castellanii*, the amoeba that Pandoravirus infects). CAI = 1 means every codon is the host-preferred choice; CAI = 0 means perfect anti-adaptation.

**RSCU distance.** Relative Synonymous Codon Usage distance — a measure of how far a gene's codon choices are from the genome-average preferred pattern. A lower RSCU distance means more genome-adapted codon usage.

**The neutrality plot.** A standard diagnostic: regress GC3 (wobble positions, neutral) against GC12 (positions 1 and 2, under selection). Under strict neutrality (no selection on codon choice), the slope should be 1.0. In *P. salinus*, the slope is β = 0.118 (95% CI: 0.089–0.146, p = 5.9 × 10⁻¹⁶). This is far below 1.0, indicating strong translational selection — codon positions under selection do not respond to changes in neutral composition the way they would if evolution were neutral. The R² of 0.045 is low, which is actually biologically informative: it means GC12 barely changes at all with GC3, consistent with very strong purifying selection at constrained positions.

**Key finding.** ORFan genes (those with no external homologs) are compositionally *indistinguishable* from annotated genes by every codon metric: GC3 difference = +0.015, CAI difference = +0.008, upstream 6-mer cosine similarity = 0.9824. This rules out recent horizontal gene transfer from external sources — a foreign gene would have compositional signatures matching its source organism, not the Pandoravirus genome.

---

## Module 3: The five-tier proto-gene continuum — an ordinal classification of gene maturity

This is the conceptual centerpiece of the paper.

**The classification scheme.** All predicted ORFs were sorted into five tiers based on two independent evidence sources: (1) whether the GenBank annotation database recognizes the sequence as a gene, and (2) whether the Prodigal *ab initio* gene predictor independently identifies it. Think of this as a 2×2 table that generates five groups based on the strength and independence of the evidence:

| Tier | Definition | n | Mean length | Mean GC |
|------|-----------|---|-------------|---------|
| T1 | Annotated + Prodigal confirmed | 516 | 1,452 bp | 0.650 |
| T2 | Hypothetical (no function) + Prodigal confirmed | 784 | 1,082 bp | 0.648 |
| T3 | GenBank-only (Prodigal missed it) | 130 | 555 bp | 0.615 |
| T4 | Prodigal high-confidence only | 335 | 501 bp | 0.560 |
| T5 | Prodigal low-confidence only | 527 | 216 bp | 0.547 |

The intergenic background (not yet a gene) sits at GC = 0.539. T5 proto-genes are at GC = 0.547 — barely above the non-coding baseline.

**The gradient.** GC content, gene length, GC3, and CAI all decline monotonically from T1 to T5. The Spearman ρ = −1.0 for both GC and length across the five tiers (exact one-tailed p = 0.0083 for each; Fisher combined p = 0.00073). A perfect rank correlation across five ordered groups is reassuring, but the biological significance is in the *consistency across independent measurements* — four different variables, measured by different methods, all agree on the same ordering. The probability that four independent measurements all align by chance is far lower than any single correlation alone.

**The proto-gene codon signal.** The most important finding in this section: even T4 and T5 proto-genes show codon optimization significantly beyond random. Compared against length-matched random intergenic ORFs:
- T4: RSCU distance 0.74 vs random 0.76 (6.7% displacement toward genome-optimized baseline; p = 0.020)
- T5: RSCU distance 1.01 vs random 1.03 (3.5% displacement; p = 0.0005)

This is important for two reasons. First, it means these proto-genes are not just random noise picked up by the gene predictor — they have genuine sequence structure. Second, because Prodigal's scoring algorithm uses hexamer (6-mer) log-likelihood at the gene start site rather than whole-gene codon optimization, the RSCU signal is partially independent of why Prodigal called them genes. The codon optimization must come from the sequences themselves.

**Spatial distribution.** Proto-genes showed zero spatial clustering near established genes (T4 p = 0.99; T5 p = 1.0 by nearest-neighbor randomization tests). This rules out mechanisms like "regulatory spillover from adjacent genes" and supports the idea that new genes can arise anywhere in the genome.

---

## Module 4: Protein structure prediction — structural relationships invisible to sequence comparison

**ESM-2 embeddings.** ESM-2 is a transformer-based protein language model (analogous to BERT for biological sequences) that generates a numerical embedding vector for any protein sequence. Two proteins with similar 3D structures tend to have similar embeddings, even if their sequences have diverged beyond the detection limit of sequence alignment tools.

We embedded 170 representative proteins spanning all five tiers and performed hierarchical clustering in embedding space, yielding 9 structural families.

**Null calibration.** A critical step: we generated composition-matched random protein sequences (preserving amino acid frequencies but randomizing order) and embedded those too. Within-family similarity among real proteins: 0.533. Among random proteins: 0.368. Cohen's d = 1.45 — a large effect size confirming the structural families are real and not artifacts of amino acid composition.

**Key finding.** 23 ORFan singletons — proteins with no detectable sequence homologs anywhere — clustered structurally with annotated domain-containing proteins. Their sequences have diverged beyond sequence-level detection, but their 3D structural organization is still recognizable. This is evidence for deep evolutionary relationships invisible to BLASTp.

**Caveat.** 41% of the variance on the first principal component correlates with GC content, indicating a partial compositional confound in the embeddings. We acknowledge this but argue the null calibration (which uses composition-matched controls) accounts for it.

---

## Module 5: The distributed AT regulatory system — the paper's central discovery

**Background on regulatory signals.** In biology, genes are not automatically expressed — the cell needs to be told where each gene starts and ends. Regulatory sequences called *promoters* (at gene starts) and *terminators* (at gene ends) serve this function. In most viruses, regulatory signals are concentrated near existing genes. The question for de novo gene birth is: how can a new gene in the middle of an intergenic region be expressed if it is far from any existing regulatory element?

**The signal.** We computed AT content across aligned gene boundaries (a "meta-gene" profile) — averaging AT content at each position relative to the gene start and end across all 1,430 genes. The signal is striking: AT content spikes by **+10.3 percentage points** at gene starts and **+9.8 percentage points** at gene ends, relative to the flanking baseline. This is happening in a genome that is 62% GC — these AT spikes are fighting uphill against the overall base composition. By comparison, Lambda shows ±1 pp and Mimivirus shows ±2 pp.

**Poly-A/T run census.** We counted all homopolymer runs of ≥ 5 bp (AAAAA, TTTTT, etc.) across the genome: 6,955 poly-A runs and 7,010 poly-T runs. Of these, **82.8–83.9% were located in intergenic regions** — far above Lambda (12–24%) and Mimivirus (16–17%). These runs are clearly concentrated outside genes.

**Strand-asymmetric positioning.** This is the key regulatory evidence. In molecular biology, promoters work on specific strands (poly-A runs upstream of genes on the coding strand; poly-T runs at the gene end on the coding strand function as terminators). We applied a strand-aware classification: for each run, we asked whether its position and strand orientation were consistent with promoter function, terminator function, or neither. Among boundary-proximal runs, **74% occupied the correct regulatory position** (binomial p < 10⁻⁶ vs. the null that correct and incorrect orientations are equally likely). This is transcriptomic validation: these AT runs are behaving like real regulatory elements.

**The critical finding: regulatory signals are uniformly distributed.** We compared AT runs that are close to gene boundaries ("boundary-proximal," n = 4,331) against AT runs dispersed throughout intergenic space (n = 7,309). Mean AT content difference: 1.6 percentage points. P-value: 0.105. Cohen's d < 0.05.

Statistically, these two populations are essentially identical. This means that the regulatory potential is not concentrated at existing gene boundaries — it is distributed uniformly throughout the entire intergenic landscape. Any piece of intergenic DNA sits near a potential promoter and terminator. A randomly emerging ORF anywhere in the genome is already flanked by the machinery it needs to be expressed.

**The proto-gene test.** T4 (upstream AT content 45.0%, p = 0.14 vs. annotated 45.9%) and T5 (46.7%, p = 0.25) proto-genes show regulatory signals statistically indistinguishable from mature annotated genes. They are born into a regulatory landscape already prepared for them.

---

## Module 6: Transcriptomic validation — checking the regulatory signal against experimental data

This module checks the regulatory signal against real experimental data, not just computational predictions.

**The external data source.** Legendre et al. (2018) performed strand-specific RNA-sequencing on Pandoravirus-infected cells — directly measuring which parts of the genome are transcribed and in which direction. Their published gene annotations incorporated this experimental data.

**Validation design.** We classified all 11,640 intergenic AT runs as either "strand-consistent with a regulatory role at the nearest gene boundary" or "strand-inconsistent." The classification requires correct strand orientation for both the run type (poly-A vs. poly-T) and its position (upstream vs. downstream of the nearest gene). 74% of boundary-associated runs passed this strict strand-aware test — far above the 50% expected if orientation were random (binomial p < 10⁻⁶).

**The hairpin comparison.** Legendre et al. (2018) proposed that palindromic RNA hairpin structures (not poly-A/T runs) are the primary termination signal, claiming ~70% of transcript ends contain hairpins. We tested this: using a strict palindrome detection algorithm (stem ≥ 6 bp, exact complementarity), observed hairpin fraction 11.6% vs. null expectation 9.7% (z = 3.6, statistically significant but only 2 pp above null). The discrepancy with Legendre et al.'s 70% likely reflects their use of energy-based RNA folding (RNAfold), which detects thermodynamically favorable structures rather than exact palindromes. Conclusion: the poly-A/T system is a real and independent termination mechanism, not just a secondary correlate of hairpin structures.

---

## Module 7: BLASTp homology analysis — measuring evolutionary age through external similarity

**What BLASTp does.** Given a protein sequence, BLASTp searches a database of essentially all known proteins (~1 billion sequences in the NCBI nr database) and finds the most similar matches using a local alignment algorithm. Results are ranked by E-value — the expected number of matches this good by chance in a database this size. A match with E-value 10⁻²⁰ is extremely unlikely to be chance; a match with E-value 0.1 is probably noise.

Our threshold: E-value < 10⁻⁵. Hits only to other Pandoravirus strains are classified as "self" (the gene is a lineage-specific orphan). Hits to anything else (any other organism) are classified as "non-self."

**The annotated gene analysis.** 521 of 528 annotated genes were successfully queried (7 errors). Of 521 successful queries:
- 426 (81.8%): orphans — hits only to other Pandoravirus strains
- 95 (18.2%): non-self hits — at least one match to another organism

The 95 non-self hits show a striking taxonomic pattern: 29.5% marine algae, 25.3% other giant viruses, 13.7% Amoebozoa (the host lineage), 11.6% bacteria, 11.6% insects. The marine algae hits reflect shared ecological habitat rather than phylogenetic proximity — Pandoravirus infects marine coastal amoebae that share habitat with marine algae. The insect hits (all ankyrin and F-box domain proteins in aphid genomes) likely reflect bacterial mobile elements (Wolbachia endosymbionts) rather than direct virus-to-insect transfer.

**The annotated vs. hypothetical comparison.** This is the key Fisher's exact test in the paper. We compared non-self hit rates between annotated genes (those with recognized protein domains) and hypothetical genes (ORFs with no recognized function):

|  | Non-self hits | Orphans | Total |
|--|--|--|--|
| Annotated | 95 | 426 | 521 |
| Hypothetical | 9 | 290 | 299 |

Fisher's exact test: OR = 7.19, 95% CI: 3.5–13.6, p = 7.79 × 10⁻¹². Annotated genes are 7.2× more likely to have external relatives than hypothetical genes. This makes intuitive sense under the maturation model: annotated genes are older (have had more time to accumulate detectable homologs in related organisms), while hypothetical genes are younger (recently born, not yet diverged to recognizable relatives in other lineages).

**E-value sensitivity analysis.** A methodological concern: what if this result is sensitive to the E-value cutoff? We repeated the annotated gene analysis at 10⁻⁵, 10⁻¹⁰, 10⁻²⁰, and 10⁻⁵⁰. Orphan rates: 81.8%, 83.1%, 84.1%, 90.4%. The conclusion is robust — if anything, tighter thresholds make orphan rates *higher*, confirming that the non-self hits are not borderline noise.

**Spatial enrichment.** Non-self hits (older genes) are enriched in the "native zone" (a genomic region from 1.3–2.0 Mb that shows hallmarks of being an older, better-annotated region): 45/95 hits fell here (OR = 1.94, Fisher's exact p = 0.004). This is spatial evidence for gene aging — the oldest genes in the genome cluster in the oldest genomic region.

**The phylogenetic distance gradient.** Hypothetical genes with non-self hits connect exclusively to phylogenetically close organisms (giant viruses, Amoebozoa). Annotated genes connect to a much broader range (marine algae, insects, diverse bacteria). Older genes have had more time to diverge and spread connections across more distantly related organisms — this is an additional prediction of the maturation model that the data confirm.

**P. dulcis exhaustive BLASTp.** We ran the same analysis on all 319 valid annotated + 741 valid hypothetical proteins of *P. dulcis* (the sister species), an exhaustive run with no subsampling. Results: 94.7% annotated orphan rate and 98.4% hypothetical orphan rate. The annotated vs. hypothetical contrast holds (OR = 3.42, p = 1.51 × 10⁻³). The *P. dulcis* annotated orphan rate (94.7%) is significantly *higher* than *P. salinus* (81.8%; OR = 0.253, p = 2.86 × 10⁻⁸), suggesting *P. dulcis* has either a younger or more lineage-isolated gene complement. The hypothetical orphan rates (98.4% vs. 97.0%) are indistinguishable — the upper bound on orphan fraction is a conserved genus-level property.

---

## Module 8: Gene family detection — adaptive GC-normalized threshold

**The problem with fixed thresholds.** A standard method for identifying duplicate genes is to compute pairwise protein sequence similarity and cluster genes above some threshold. In a 62% GC genome, protein amino acid composition is skewed toward G/C-coded amino acids (Ala, Gly, Pro, Arg), which inflates pairwise similarity among *unrelated* proteins just because they have similar letter frequencies. In an AT-rich genome like Mimivirus (28% GC), the bias runs the other direction.

Using a fixed threshold of cosine similarity > 0.95 on 3-mer amino acid frequencies: Mimivirus null distribution has a mean of 0.9593, meaning *random unrelated proteins* score above 0.95. The fixed threshold classifies essentially all Mimivirus protein pairs as similar — a catastrophic false positive rate.

**The adaptive threshold.** For each genome independently, we generated 5,000 pairs of random sequences with the genome's actual base composition and gene length distribution, computed pairwise 3-mer cosine similarities, and set the clustering threshold at the 99th percentile of this null distribution. This ensures a consistent ~1% false-positive rate across all genomes regardless of base composition.

**Results.** Under adaptive thresholds:
- *P. salinus*: 64.1% singletons (no similar partner found)
- *P. dulcis*: 68.7%
- Mimivirus: 88.9% (vs. the artifact-driven 9.5% under the fixed threshold)
- Pithovirus: 90.4%
- Mollivirus: 81.6%

All five giant virus genomes are singleton-dominated — gene duplication is not the primary engine of genome expansion. Duplication is a minor contributor; de novo origin is the dominant mode.

**Sensitivity analysis.** Repeating at the 95th and 99.5th percentile thresholds rather than the 99th: singleton rates shift by ~2–4 percentage points across all genomes. The qualitative conclusion — singleton-dominated proteomes — is robust.

---

## Module 9: Cross-NCLDV comparative analysis

**Design.** We applied Modules 2–7 to all five NCLDV genomes (spanning 28–64% GC, 610 kb to 2.5 Mb, four virus families) and looked for which patterns are consistently conserved vs. which are lineage-specific.

**Three universally conserved patterns:**

1. **Strand balance near 50%.** Four of five genomes have 48.8–51.3% of genes on the plus strand. *Pithovirus* is the outlier at 44.3% (approximately 2.5 standard deviations from 50% given its gene count), possibly reflecting replication-associated strand asymmetry. The near-perfect strand balance in the other four is consistent with no strong directional pressure on gene orientation.

2. **Translational selection toward host codon usage.** Neutrality plot slopes range from 0.067 to 0.173 (all well below the 1.0 expected under neutrality) for four of five genomes, indicating selection is pushing codons toward host-optimized usage in all cases. *Mollivirus* is an anomaly with a negative slope (−0.146), meaning GC3 and GC12 move in opposite directions — an unusual pattern whose cause is unclear (possibly annotation quality or genuine unusual selection).

3. **Singleton-dominated proteomes.** Under adaptive thresholds, 64–91% singletons across all five genomes. De novo gene birth, not duplication, is the dominant growth strategy.

**The maturation direction reversal — the strongest cross-genome finding.** In all four GC-rich genomes (*P. salinus*, *P. dulcis*, Mollivirus, Pithovirus), annotated genes have *higher* GC3 and CAI than hypothetical genes — maturation is driving codon usage toward the high-GC genome preference. In Mimivirus (28% GC), the direction *reverses*: annotated genes have *lower* GC3 (Δ = −0.024, Mann-Whitney p = 1.09 × 10⁻⁹, Cohen's d = −0.439) and lower CAI (Δ = −0.012, p = 2.30 × 10⁻¹³, d = −0.467) than hypothetical genes. This is exactly what the maturation model predicts: older genes are more adapted to *host* codon preferences. In a high-GC host, host-adapted = more GC. In a low-GC host, host-adapted = more AT. The reversal is a genome-composition-independent confirmation of the model.

If maturation didn't drive toward host adaptation, if instead genes randomly drifted compositionally, you would not expect systematic reversal based on genome GC content. The reversal is a falsifiable, verified prediction.

---

## The four-stage maturation model — integrating everything

The paper's synthesis of all findings is a proposed temporal ordering for de novo gene birth:

**Stage 1: Regulatory signals (pre-existing).** Poly-A/T regulatory signals pre-exist throughout intergenic space (82.8–83.9% of runs are intergenic; boundary-proximal and dispersed runs are compositionally identical; proto-genes T4/T5 already have mature regulatory signatures). No action needed — any stretch of intergenic DNA already has promoter/terminator infrastructure.

**Stage 2: Protein structure.** Proto-genes cluster with established protein structural families in embedding space (nearest-neighbor similarity test p = 0.90, meaning proto-genes are not structurally distinguishable from established genes). Structural organization emerges early, even before sequence-level recognition.

**Stage 3: Codon optimization.** T4 and T5 proto-genes show statistically significant displacement toward genome-preferred codon usage relative to random ORFs (p = 0.020 and p = 0.0005). This signal is detectable in the least mature candidates, indicating selection has begun to optimize codon choices.

**Stage 4: GC composition.** This is the slowest process. Proto-genes (GC = 0.547) have barely moved from the intergenic baseline (GC = 0.539) toward the coding gene average (GC = 0.650). Full compositional equilibration is the last step and likely takes the most evolutionary time.

**Critical caveat.** This ordering is inferred from a **cross-sectional snapshot** — we are looking at genes at different maturation stages simultaneously, not tracking individual genes through time. This is exactly analogous to cross-sectional epidemiology vs. longitudinal studies. The cross-sectional gradient is consistent with the proposed ordering but cannot prove necessity (i.e., that Stage 2 cannot precede Stage 1 for any individual gene). A definitive test requires comparative analysis across multiple Pandoravirus strains at different evolutionary distances from each other — the longitudinal equivalent.

---

## What this paper is *not* claiming

**It is not claiming proof of mechanism.** The paper provides computational evidence consistent with a model. The four-stage ordering is a hypothesis, not a demonstrated causal chain.

**It is not claiming the regulatory signals are proven promoters/terminators.** The 74% strand-consistent positioning and the +10.3 pp AT spikes are strongly suggestive, but direct experimental confirmation (reporter assays, TSS/TTS mapping) has not been done. The transcriptomic validation uses annotation-informed GFF files, not raw read-level TSS coordinates.

**The ESM-2 structural clustering is a proxy.** The model generates structural embeddings, not actual 3D protein structures. It is consistent evidence, not structural proof.

**Multiple testing is acknowledged.** 12 analytical modules, many statistical tests. The paper distinguishes between four pre-specified confirmatory tests (the ones that directly test the central hypothesis) and the many exploratory comparisons. The four confirmatory tests are: (1) proto-gene codon usage vs. random ORFs, (2) the strand-aware 74% binomial test, (3) the annotated-vs-hypothetical Fisher's exact test, and (4) the two-population AT run test.

---

## Why this is interesting to a statistician

Beyond the biology, this paper has several methodologically notable features:

1. **Effect-size discipline in genomics.** The explicit rejection of 88.6% of nominally significant k-mer signals based on effect-size filters is unusual in the genomics literature, where z > 3 is often treated as sufficient. The paper makes the case that in large genomes, statistical power is not a meaningful filter.

2. **Composition-aware null models.** The adaptive GC-normalized threshold for gene family detection solves a real confounding problem that invalidates a widely-used methodology in AT-rich genomes. The old Mimivirus singleton rate of 9.5% (a false-positive artifact) vs. the corrected 88.9% shows the stakes of ignoring this.

3. **Cross-sectional vs. longitudinal inference.** The explicit acknowledgment that the maturation ordering is cross-sectional, and the statement that longitudinal confirmation requires comparative multi-strain tracking, reflects a statistical rigor that is often glossed over in biological papers that present cross-sectional gradients as temporal sequences.

4. **Distinguishing confirmatory from exploratory tests.** The paper explicitly separates its four pre-specified confirmatory tests from the broader exploratory analysis. This is methodological transparency that is more common in clinical trials literature than genomics.

5. **The biological "reversal" as a falsifiable prediction.** The codon maturation direction reversal in Mimivirus is the kind of result that is hard to obtain by chance or confounding — it requires the direction of a biological process to switch based on an independent variable (genome GC content). This is close to a natural experiment, and the result is exactly what the model predicts.

---

*Document prepared March 2026. All cited statistics are mapped to source computations in statistics_audit.json.*
