**Eddie Bradford**
Independent Researcher
Eddie.Bradford@gmail.com
ORCID: 0009-0005-5870-1066

March 28, 2026

Editorial Board
*Genome Biology and Evolution*
Oxford University Press

---

Dear Editors,

I am writing to submit the manuscript **"Computational evidence for a distributed AT-based regulatory mechanism enabling de novo gene birth in *Pandoravirus salinus*"** for consideration as a Research Article in *Genome Biology and Evolution*.

## The question

*Pandoravirus salinus* carries the largest known viral genome at 2.47 Mb, with an estimated 91.4% (95% CI: 89.6%–93.0%) of its ~1,430 predicted genes classified as ORFans — sequences with no detectable homologs outside the Pandoraviridae. While the Claverie group has established through comparative genomics that de novo gene birth from non-coding sequence is the most likely explanation, the genomic infrastructure enabling this process at scale has remained uncharacterized. This manuscript provides the first systematic, multi-module computational analysis of that infrastructure in a single genome, with comparative validation across four additional Nucleocytoviricota spanning four families and a 36-percentage-point range in GC content.

## Key findings

**A pre-existing regulatory scaffold drives the core argument.** Using strand-aware analysis, we show that poly-A/T runs with the positional and orientation hallmarks of eukaryotic transcriptional regulatory signals are present throughout *P. salinus* intergenic space at gene boundaries (AT enrichment +10.3/+9.8 pp upstream/downstream; 74% of boundary-associated runs in correct regulatory position, p < 10⁻⁶). These signals pre-exist gene formation, providing the mechanistic substrate for transcription of nascent ORFs regardless of genomic position. Comparative analysis of *Pandoravirus dulcis* confirms this distributed AT regulatory system operates at the genus level: 86% intergenic poly-A/T fraction, +7.7 pp gene boundary AT enrichment, and 61.8% correct strand-aware regulatory positioning — fully consistent with *P. salinus* values.

**A quantitative five-tier proto-gene continuum.** We identify a continuous gradient from intergenic sequence through four hypothetical-gene classes to established annotated genes, with monotonic trends across GC content, length, and codon adaptation index. ESM-2 protein language model embeddings reveal nine structural families that are cross-tier — structural relationships invisible to sequence homology methods — consistent with early functional divergence preceding sequence consolidation. *P. dulcis* replicates the same five-tier compositional architecture independently, transforming this from a single-species observation to a Pandoraviridae-level property.

**Two null calibrations with opposed effect sizes.** A key methodological strength is the direct quantification of two competing artifact explanations:

- *Gene duplication artifact*: Gene family analysis shows the diffuse pairwise similarity cloud (mean 0.109 vs Mimivirus 0.069) is fully explained by *P. salinus*'s GC-skewed amino acid composition. A composition-matched null model reproduces the background (null mean 0.105; **Cohen's d = 0.06** — near-zero effect). The apparent similarity is compositional noise, not duplicated sequence.
- *Structural signal artifact*: ESM-2 structural families are genuine. Within-family embedding similarity (0.533) far exceeds the composition-matched null (0.368; **Cohen's d = 1.45** — very large effect). Structural signal is real where sequence signal is not.

The contrast between these two calibrations — d = 0.06 vs d = 1.45 within the same dataset — provides strong internal validation that the analysis discriminates artifact from signal.

**A four-stage maturation model confirmed across five NCLDVs.** Converging evidence supports the sequence: regulatory signals (pre-existing) → nascent protein structure → codon optimization → GC equilibration. Exhaustive BLASTp analysis of all *P. salinus* annotated genes (n = 521) and a large hypothetical gene subsample (n = 299) shows annotated genes are substantially more likely than hypothetical genes to have non-Pandoravirus homologs (OR = 7.19, 95% CI: 3.5–13.6, p = 7.79 × 10⁻¹²), providing an age-dependent homology gradient consistent with the maturation model. Formal Mann-Whitney U tests for GC3 and CAI maturation shifts confirm the model in all five NCLDV genomes tested (all p < 0.05), with the direction of codon maturation reversing in AT-rich Mimivirus (28% GC) relative to the four GC-rich genomes — precisely as the de novo birth model predicts. Exhaustive BLASTp analysis of *P. dulcis* (all 323 annotated + all 747 hypothetical genes) confirms genus-level consistency: annotated orphan rate 94.7%, hypothetical orphan rate ~98%, with a Fisher's exact test comparing *P. dulcis* vs *P. salinus* annotated orphan rates yielding a statistically significant difference (p < 0.001), suggesting *P. dulcis* annotated genes are modestly more lineage-specific than their *P. salinus* counterparts — an interesting evolutionary divergence within the genus.

**Cross-NCLDV comparative validation.** Three features are consistent across all five NCLDV genomes tested: near-perfect strand balance (48.8–51.3% plus-strand in four of five genomes), translational selection toward host codon usage, and singleton-dominated proteomes (64–91% singletons with adaptive thresholds). The AGCT tetranucleotide ranks as the single most depleted 4-mer of all 256 (rank 1/256) in both Pandoravirus species and in *Mollivirus sibericum*, while enriched in *Mimivirus* (rank 241/256) and *Pithovirus* (rank 185/256) — implicating a lineage-specific sequence constraint in the Pandoraviridae beyond what GC content alone predicts.

## Methodological transparency

This manuscript was developed under an explicit pre-submission statistical audit protocol. Several specific transparency measures are worth noting to reviewers:

1. **ORFan rate derived independently.** Rather than citing Philippe et al. (2013)'s 92.5% as our own result, we compute our estimate directly from BLASTp data using an exhaustive annotated sample and large hypothetical subsample. Philippe et al.'s figure serves as independent confirmation by a complementary method.

2. **Adaptive GC-normalized gene family detection.** Fixed cosine similarity thresholds fail for AT-rich genomes (Mimivirus null mean 0.9593 exceeds the 0.95 fixed threshold). We developed an adaptive threshold method using genome-matched random sequence null distributions (99th percentile), enabling valid cross-genome comparison. The old Mimivirus singleton rate of 9.5% was a method artifact; the corrected rate (88.9%) makes the five-genome picture coherent.

3. **Haldane-corrected odds ratios throughout.** All 2×2 contingency statistics report both Fisher's exact OR and the Haldane-corrected OR with 95% CI.

4. **Dinucleotide null validation.** Markov-1 dinucleotide-preserving null models produce gene family thresholds differing by only Δ = +0.0001 from mononucleotide null models, confirming that mononucleotide null calibration is adequate.

5. **Machine-readable audit trail.** A `statistics_audit.json` file included in the repository maps every manuscript-cited number to its source computation, including confidence intervals, sample sizes, and test specifications.

All analysis scripts, the PostgreSQL database schema, intermediate results, BLAST output, and the `statistics_audit.json` are available at https://github.com/mathamagician/pandoravirus-denovo-genebirth and archived with a permanent DOI at https://doi.org/10.5281/zenodo.19046142.

## Fit to GBE

*Genome Biology and Evolution* is the natural venue for this work: it addresses a central question in genome evolution (de novo gene origin) through rigorous computational analysis of a landmark genome, with comparative validation across four additional giant virus families. The subject organism was first described in *Science* and has attracted ongoing interest from the evolutionary genomics community. The 12-module framework, adaptive null model methodology, and transparent audit approach may also be of methodological interest to readers working on other orphan-gene-rich genomes.

This manuscript has not been published and is not under consideration elsewhere. AI tools (Anthropic's Claude) were used to assist with code development and editorial refinement, as disclosed in the manuscript. All analyses, statistical methods, interpretations, and final content were reviewed and validated by the author.

I look forward to your consideration.

Sincerely,

**Eddie Bradford**
Independent Researcher
Eddie.Bradford@gmail.com
ORCID: 0009-0005-5870-1066

---

*Suggested reviewers:*

- **Jean-Michel Claverie** (Aix-Marseille University) — original Pandoravirus discoverer; deep expertise in Pandoraviridae genomics
- **Natalya Yutin** (NCBI) — giant virus comparative genomics, de novo gene evolution
- **Eugene Koonin** (NCBI) — viral gene origin, de novo gene birth theory
- **Sébastien Santini** (CNRS/IGS Marseille) — bioinformatics co-author on core Pandoravirus genome papers
- **Betül Kaçar** (University of Wisconsin) — molecular evolution, ancestral sequence reconstruction

*Reviewers to exclude:* None requested.
