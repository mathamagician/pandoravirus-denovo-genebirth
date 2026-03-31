# Statistical Methodology Audit — Pandoravirus 1.0 + 2.0 Combined Paper
**Document:** `Pandoravirus_combo_paper_draft3_GBE.md`
**Audit date:** 2026-03-21
**Auditor:** Principal Statistician perspective (RAM/RPSM/consult-stats framework)
**Target journal:** Genome Biology and Evolution (Oxford)
**Purpose:** Pre-submission statistical audit — identify errors, weaknesses, and improvements before peer review

---

## Executive Summary

The paper is statistically stronger than average for a single-investigator computational genomics study. The use of multiple null calibrations (Phases 7b and 8b), effect-size filtering alongside p-values, Fisher's exact with Haldane-corrected Woolf CIs, and an explicit 140-entry statistics audit JSON are genuine differentiators. The peer review process already caught and corrected two real errors (taxonomy misclassification, connection-failure/no-hit conflation).

**RPSM-style pre-submission score: ~64/100 (Provisional — ⭐⭐⭐)**
- S_stat: 62/100 — several tests under-specified; multiple testing unaddressed; key numbers without CIs
- S_evid: 55/100 — single genome, N=4 comparison ladder, no wet lab validation
- S_transp: 90/100 — GitHub repo, Zenodo DOI, statistics_audit.json, code available
- S_rep: 50/100 — new paper, unverified; partially consistent with Claverie group findings

This is a **publishable paper** at GBE with revisions. The issues below are the ones a GBE reviewer is most likely to catch. Fix Tier 1 before submission; address Tier 2 in the manuscript to pre-empt reviewer objections; Tier 3 items strengthen credibility but are not submission blockers.

---

## TIER 1 — Must Fix Before Submission

These issues could generate major revisions or rejection if found by a reviewer.

---

### 1.1 The 92.5% ORFan rate is not derivable from the reported data

**Where:** Abstract ("92.5% of genes classified as ORFans"), Limitations section.

**The problem:** From the raw BLASTp numbers reported in the paper:
- Annotated orphans: 426 / 521 queried = 81.8%
- Hypothetical orphans: 290 / 299 queried = 97.0%

Extrapolating to all 1,430 genes:
- Annotated estimate: 81.8% × 528 = 432
- Hypothetical estimate: 97.0% × 902 = 875
- **Total: 1,307 / 1,430 = 91.4%** — not 92.5%

If all 902 hypothetical genes are counted as orphans (100% assumption): 432 + 902 = 1,334 / 1,430 = **93.3%** — also not 92.5%.

The number 92.5% does not arise from any defensible combination of the reported statistics. It appears to be a rounding estimate or was computed from an earlier version of the data.

**Fix (two options):**
- **Option A (preferred):** Report 91.4% with explicit derivation: "BLASTp-estimated orphan fraction of all 1,430 genes, applying observed rates of 81.8% (annotated) and 97.0% (hypothetical) to the full gene counts." Also provide a bootstrap CI on this extrapolation.
- **Option B:** Report 87.3% (the directly observed rate from queried genes: 716/820) and note that this is a lower bound due to the random hypothetical subsample.

Either option requires updating the abstract, the limitations section, and any other occurrence of "92.5%." This is a headline number and should be exact and derivable.

**RESOLUTION (draft4, 2026-03-21):** Fixed. Abstract updated to: "estimated 91.4% (95% CI: 89.6%–93.0%) of genes classified as ORFans, based on our BLASTp analysis (annotated-gene orphan rate 81.8%; hypothetical-gene orphan rate 97.0%; extrapolated to all 1,423 classified genes; Philippe et al. 2013 independently reported 92.5%)." Introduction updated to cite Philippe 2013 for 92.5% and note our independent 91.4% estimate. Limitations updated. statistics_audit.json updated with derivation. Bootstrap CI computed via compute_missing_statistics.py.

---

### 1.2 Neutrality plot: slope not reported, only R²

**Where:** Module 4, Results section; Methods.

**The problem:** The neutrality plot (GC12 vs GC3) is interpreted as showing "translational selection dominant, not mutational drift." This interpretation is based on **R² < 0.05** (weak correlation). But R² alone does not tell you what you need to know. The correct interpretation depends on the **slope** of the GC12-on-GC3 regression:
- Slope ≈ 0: GC12 is invariant while GC3 varies → selection fixes constrained positions regardless of mutational pressure → selection dominant ✓
- Slope ≈ 1: Both GC12 and GC3 vary together → mutation dominant
- R² = 0.05 with any slope: noisy relationship; it's the slope that determines direction

R² < 0.05 is consistent with slope ≈ 0 (selection) but also with a noisy slope > 0. A GBE reviewer will ask for the slope.

**Fix:** Report the linear regression slope (GC12 ~ GC3), 95% CI on the slope, and the p-value. The interpretation "selection dominant" only follows from a slope near zero with a narrow CI. Add one sentence to Results: "The regression slope (β = X, 95% CI: Y–Z) indicates that GC content at constrained positions is largely invariant with changes in wobble-position GC, consistent with selection rather than mutational drift."

**RESOLUTION (draft4, 2026-03-21):** Fixed. Computed via gene_codon table: slope = 0.1176 (95% CI: 0.0894–0.1457), R² = 0.0448, p = 5.9×10⁻¹⁶, N = 1,430 genes. Added to Results codon usage paragraph. Slope << 1.0 confirms selection dominant over mutational drift.

---

### 1.3 CAI comparison lacks any statistical test or CI

**Where:** Module 4, Results; Abstract.

**The problem:** "Pandoravirus CAI 0.787 vs Mimivirus 0.657 vs host 1.000" is stated as if these are exact population values. In reality, these are means over 1,430 and 979 per-gene CAI distributions. The comparison ("Pandoravirus is better host-adapted than Mimivirus") is a central claim in the paper and requires a formal test.

With 1,430 Pandoravirus genes and 979 Mimivirus genes, a Mann-Whitney U test is trivial to run and will almost certainly give p << 0.001. More importantly, **report the per-gene standard deviations and a 95% CI on the difference in means** — this lets the reader evaluate the magnitude of the difference.

**Fix:** Add to Results: "Mean per-gene CAI was 0.787 ± [SD] for Pandoravirus vs 0.657 ± [SD] for Mimivirus (Mann-Whitney U, p = X; difference = 0.130, 95% bootstrap CI: Y–Z)."

**Secondary issue:** The CAI caveat (Kazusa general codon table, not highly-expressed gene subset) is stated in Limitations but should also appear at first mention in Results. GBE reviewers are sophisticated — they'll catch this immediately and will want to see it proactively disclosed at the point of use.

**RESOLUTION (draft4, 2026-03-21):** Fixed. Results updated with full statistics: P. salinus mean CAI 0.787 ± 0.032 (N=1,430) vs Mimivirus 0.657 ± 0.027 (N=979); Mann-Whitney U = 1,392,180, p < 10⁻³⁰⁰; mean difference 0.129, 95% bootstrap CI: 0.127–0.132. CAI reference caveat embedded at point of use.

---

### 1.4 OR = 11.9 for hypothetical spatial enrichment lacks a CI

**Where:** Module 11-12, Results (Homology landscape section).

**The problem:** "7/9 hypothetical hits in native zone (OR = 11.9, p = 0.0009)" is reported without a CI. With n = 9 total hypothetical hits (7 in native zone, 2 outside), Fisher's exact test with Haldane correction should produce an extremely wide CI — likely something like OR = 11.9, 95% CI: 1.8–1,500 or wider. The point estimate of OR = 11.9 looks impressive but is based on a 2×2 table with very small cell counts.

Omitting the CI when n is this small is a presentation problem — a reviewer will notice immediately that OR = 11.9 from n = 9 is not robust. Reporting the CI pre-empts this objection and makes the finding more credible by showing honest uncertainty.

**Fix:** Calculate and report the Fisher's exact 95% CI for OR = 11.9. If the CI lower bound is reasonably above 1.0 (say > 1.5), the finding is still presentable. If the CI is extremely wide, frame it as "preliminary evidence" or move to supplementary. Either way, the CI must be reported.

**RESOLUTION (draft4, 2026-03-21):** Fixed. Computed contingency table: hits_in_nz=7, nohits_in_nz=67, hits_out_nz=2, nohits_out_nz=224. Fisher OR = 11.70, p = 0.0010. Haldane-corrected OR = 9.98, 95% CI: 2.32–42.83. CI is wide as expected (n=9 hits) but lower bound 2.32 > 1.0. Added to Results with framing: "the wide CI reflects small n=9, but even the lower bound 2.32 indicates substantial enrichment." Note: the paper originally cited OR=11.9, which was a prior rough calculation; corrected value is 11.70.

---

### 1.5 Multiple testing across 12 modules never explicitly addressed

**Where:** Methods, Limitations.

**The problem:** The paper reports ~50+ statistical tests across 12 modules: k-mer z-scores (thousands of k-mers), per-gene codon comparisons, several Fisher's exact tests, binomial tests, Mann-Whitney tests, Spearman correlations, chi-squared tests. No Bonferroni, Holm, or Benjamini-Hochberg correction is ever mentioned.

GBE reviewers will ask: "With this many tests, how many findings would you expect by chance?" The current paper implies every p < 0.05 result is meaningful, but with ~50 tests and no correction, you'd expect ~2.5 false positives by chance even if all nulls were true.

**The correct response is NOT to apply Bonferroni across all tests** (they test very different hypotheses and are not an experiment-wise family). The correct response is to:
1. Explicitly acknowledge the multiple testing situation in the Limitations or Methods
2. Clarify which tests are **confirmatory** (pre-specified, part of the central hypothesis) vs **exploratory** (suggested by the data)
3. Note that the central claims are supported by **convergent evidence from independent methods** — this is a genuine protection against multiple testing

**Fix:** Add one paragraph to Limitations: "Multiple statistical tests were performed across 12 analytical modules. We did not apply a global multiple testing correction because the tests address distinct biological questions at different scales; application of Bonferroni correction across unrelated tests would be overly conservative and inappropriate. Instead, we distinguish between confirmatory tests (the proto-gene codon usage test, the strand-aware binomial test, and the annotated-vs-hypothetical Fisher's exact test, each pre-specified before data analysis) and exploratory comparisons that generated the hypotheses tested by the confirmatory analyses. The central conclusions — that P. salinus has a distributed regulatory system and a measurable maturation gradient — are supported by multiple independent lines of evidence rather than by any single p-value."

**RESOLUTION (draft4, 2026-03-21):** Fixed. Added full multiple testing paragraph to Limitations, distinguishing confirmatory tests (proto-gene codon usage Mann-Whitney, strand-aware binomial, annotated-vs-hypothetical Fisher's) from exploratory comparisons. Added within-genome independence caveat. Explicitly states conclusions rest on convergent large-effect-size evidence.

---

### 1.6 Prodigal overlap criterion is ambiguously defined

**Where:** Module 5 Methods; Results.

**The problem:** The methods state "≥80% reciprocal overlap matching" but does not specify:
- 80% of which gene's length? (the longer one? shorter one? both?)
- Does start position need to match, or just overlap?

This matters for the novel prediction count (347 Prodigal-only predictions) and for the "true intergenic drops from 32.5% to 23.8%" calculation. Different definitions can produce substantially different counts.

**Fix:** Specify precisely: "A Prodigal prediction was considered to confirm a GenBank annotation if at least 80% of the GenBank annotation's length overlapped with the Prodigal prediction and at least 80% of the Prodigal prediction's length overlapped with the GenBank annotation (reciprocal 80% overlap). Start position matching was not required."

**RESOLUTION (draft4, 2026-03-21):** Fixed. Methods updated to: "A Prodigal prediction was considered to confirm a GenBank annotation if (overlap length / GenBank annotation length) >= 0.80 AND (overlap length / Prodigal prediction length) >= 0.80; both conditions must be met simultaneously, and start-position matching was not required."

---

### 1.7 Proto-gene codon usage test: test statistic never specified

**Where:** Module 7c Methods and Results.

**The problem:** "T4 codon usage closer to genome than random IG ORFs (p = 0.020, T5 p = 0.0005)" — but what test? Mann-Whitney? Permutation test? KS test? On what distance metric exactly? The RSCU Euclidean distances are mentioned, but the test that generated p = 0.020 and p = 0.0005 is never identified in the methods section.

GBE requires that methods be reproducible. A reader cannot reproduce this result without knowing the test.

**Fix:** Add to Methods (Module 7c): "For each proto-gene in T4 and T5, RSCU Euclidean distance to the genome-preferred codon vector (mean RSCU across all T1+T2 genes) was computed. The same distance was computed for 10,000 length-matched random intergenic ORFs. Mann-Whitney U test compared observed T4/T5 distances against the random null distribution."

**RESOLUTION (draft4, 2026-03-21):** Fixed. Methods (Module 7c) updated with full specification: RSCU Euclidean distance to T1+T2 mean vector, 10,000 length-matched random intergenic ORFs null, one-sided Mann-Whitney U test.

---

## TIER 2 — Should Address (Pre-empts Reviewer Objections)

These issues don't necessarily block publication but will attract reviewer questions. Addressing them proactively signals statistical sophistication.

---

### 2.1 N = 4 genome comparison — all "strongest" claims are limited to these 4 genomes

**Where:** Abstract; Results (throughout).

**The problem:** Claims like "the lowest GC coefficient of variation," "the strongest coding-intergenic GC gap (10.7 pp)," and "the strongest gene boundary signals tested" are presented as if Pandoravirus is uniquely extreme among viruses. In reality, these are comparisons to **3 other genomes chosen to span the dsDNA size spectrum** — PhiX174 (5 kb), Lambda (48 kb), Mimivirus (1.2 Mb). This is a non-random, deliberately informative comparison set, not a representative sample of known viruses.

A GBE reviewer could ask: "Have you compared to the other 14 known Pandoravirus strains? To other Nucleocytoviricota? To other giant DNA viruses?" Among all known Pandoravirus genomes, is P. salinus the most compositionally uniform? Among all characterized viruses with poly-A/T regulatory systems, is the +10.3/+9.8pp spike the strongest?

**Fix:** Every comparative superlative should be qualified with "among the comparison genomes tested" or "in this four-genome comparison." The Methods should explicitly explain the comparison ladder logic: "We selected three comparison genomes spanning the dsDNA virus size spectrum to contextualize P. salinus — one small phage (PhiX174, 5 kb), one medium phage (Lambda, 49 kb), and one other giant DNA virus (Mimivirus, 1.2 Mb). These comparisons are illustrative rather than statistical; characterization of the full Pandoraviridae family awaits larger-scale comparative analysis."

The "strongest" claims then become clearly bounded.

**RESOLUTION (draft4, 2026-03-21):** Fixed. Abstract updated to "among the four genomes tested." All comparative superlatives in Results qualified with "among the four comparison genomes" or "in our comparison set." Methods note added: "These comparisons are illustrative rather than exhaustive; characterization across the full Pandoraviridae family awaits larger-scale analysis."

---

### 2.2 GMM clustering on highly correlated features — GC confound risk

**Where:** Module 6 Methods; Results (evolutionary origin classification).

**The problem:** The GMM uses 5 features: GC, GC3, GC12, ENC, CAI. These 5 features are all substantially driven by a single underlying signal — genome-wide GC content (GC → GC3, GC3 influences GC12, GC composition influences amino acid composition → CAI). Running GMM on 5 correlated features that all load heavily on the same factor inflates apparent cluster separation (because the cluster separation is measured in a 5D space where 4+ dimensions carry the same signal).

Phase 8 already identified that 41% of PC1 variance in ESM-2 embeddings correlates with GC — the same issue exists in Module 6 but was never quantified there. If GC alone accounts for most of the cluster separation, the "3 clusters" are really 3 GC-content bands, not 3 evolutionary origins.

**Fix (level A — paper-level):** Add to Limitations: "The GMM clustering in Module 6 uses five features (GC, GC3, GC12, ENC, CAI) that are highly correlated with genome-wide GC content. The three clusters identified may primarily reflect quantitative GC variation rather than distinct evolutionary origins; interpretation as 'NATIVE,' 'HIGH-GC,' and 'VARIANT' populations should be considered provisional pending analysis with GC-corrected or dimensionality-reduced feature spaces."

**Fix (level B — if time permits before submission):** PCA-reduce the 5 features to 2-3 orthogonal components, re-run GMM on PCA scores, compare the resulting clusters to the current 3-cluster solution. If results differ materially, report the PCA-GMM as the primary analysis and the original as supplementary.

**RESOLUTION (draft4, 2026-03-21):** Level A fix applied. Added to Limitations: "The GMM clustering uses five features highly correlated with genome-wide GC content; the three clusters may primarily reflect quantitative GC variation; designations should be considered provisional."

---

### 2.3 Annotation-based proto-gene tiers create partial circularity

**Where:** Module 7b; Results and Discussion.

**The problem:** The five-tier continuum (T1-T5) is defined by annotation confidence:
- T1: GenBank annotated + Prodigal confirmed
- T2: GenBank hypothetical + Prodigal confirmed
- T3: GenBank-only (Prodigal missed)
- T4: Prodigal-only, high confidence
- T5: Prodigal-only, low confidence

The gradient (GC 0.650 → 0.547, length 1,452 → 216 bp) is interpreted as evidence for a maturation continuum. BUT: NCBI annotation practice assigns lower confidence to shorter, lower-GC sequences, and Prodigal preferentially calls shorter ORFs in more AT-rich regions with lower GC. So the gradient in GC and length across tiers could be **partially driven by annotation practice**, not by biological maturation. The shorter, lower-GC sequences are placed in lower tiers BECAUSE annotation tools are less confident about them — the tiers predict their own GC/length distribution.

This is an annotation feedback loop (related to the RAM discussion): the annotation tier is used as a proxy for biological maturity, but the annotation tier itself is based on sequence properties (length, GC) that create the gradient.

The paper acknowledges this ("the gradient could reflect annotation quality") but does not attempt to quantify how much of the gradient is attributable to annotation bias vs biological signal. The Phase 7c disambiguation tests (codon usage and spatial proximity) partially address this but do not close the gap.

**Fix:** Add to Discussion: "The gradient in GC and length across annotation tiers is consistent with a biological maturation model, but we cannot exclude the alternative that annotation practice contributes to the observed gradient — shorter, lower-GC sequences are inherently harder to annotate, and Prodigal may preferentially detect shorter ORFs in AT-richer intergenic regions. The Phase 7c codon usage test argues against a pure annotation artifact (T4/T5 sequences have measurably different codon preferences from length-matched random ORFs, indicating sequence-level structure), and the spatial randomness of proto-gene positions argues against operon extension. However, a definitive separation of annotation bias from biological signal requires comparative genomic validation in other Pandoravirus strains."

**RESOLUTION (draft4, 2026-03-21):** Fixed. Added circularity acknowledgment to Discussion's alternative models section: annotation practice assigns lower confidence to shorter, lower-GC sequences; Phase 7c codon usage test is strongest counter-evidence; comparative genomic validation needed for definitive separation.

---

### 2.4 The proto-gene "codon usage is 10-15% of the way from random to genome-optimized"

**Where:** Results; Discussion; Figure 1 legend.

**The problem:** This effect size description is qualitative. Looking at the RSCU distances in the findings report: Genome-optimized (T1+T2) RSCU distance = 0.46, T4 = 0.74, T5 = 1.01, Random IG ORFs ≈ 0.76 / 1.03.

T4 distance is 0.74 vs random 0.76. The "displacement" is (0.76 - 0.74) / (0.76 - 0.46) = 0.02 / 0.30 = **6.7%**, not 10-15%. T5 distance 1.01 vs random 1.03: displacement = (1.03 - 1.01) / (1.03 - 0.46) = 0.02 / 0.57 = **3.5%**, also not 10-15%.

If the 10-15% figure is calculated differently (e.g., using a different distance metric or a different baseline), the calculation should be shown. The claim "10-15% displacement" as written appears to overstate the effect size.

**Fix:** Show the explicit calculation in the paper or supplementary. If the actual displacement is ~5-7%, report that number honestly. A 5% displacement that is statistically significant (p = 0.0005) is still meaningful — codon optimization is expected to be weak in proto-genes. Overstating the effect size actually weakens the paper if a reviewer recalculates it.

**RESOLUTION (draft4, 2026-03-21):** Fixed. Confirmed from phase7c_disambiguation.py output: GenBank baseline = 0.46, T4 distance = 0.74 vs random 0.76 → displacement = (0.76-0.74)/(0.76-0.46) = 6.7%; T5 = 1.01 vs random 1.03 → displacement = 3.5%. All three occurrences of "10-15%" corrected to "3.5–6.7%" with explicit RSCU distances in text. Qualitative conclusion unchanged: effect is small but statistically significant and real.

---

### 2.5 MITE transposable elements as a source of poly-A/T runs — alternative not addressed

**Where:** Discussion (Additional compositional forces section).

**The problem:** Sun et al. (2015) showed that MITE transposable elements colonize P. salinus intergenic space. MITEs are typically AT-rich and frequently terminate with poly-A/T tails. The paper cites this work but does not address the obvious alternative explanation: **the poly-A/T runs may be residual MITE sequences rather than functional regulatory signals**. If most poly-A/T runs are in MITE footprints, the regulatory signal may be a byproduct of transposon activity rather than a designed regulatory system.

This is a potentially important alternative that a GBE reviewer will raise if the authors don't.

**Fix:** Add 2-3 sentences to the Discussion: "An important alternative source of poly-A/T runs in intergenic space is residual sequence from MITE transposable elements, which are AT-rich and common in P. salinus intergenic regions (Sun et al. 2015). We cannot currently distinguish MITE-derived poly-A/T runs from those associated with independent regulatory function using sequence analysis alone. However, several observations are consistent with an active regulatory role beyond MITE residues: the strand-asymmetric positioning of poly-A and poly-T runs at promoter and terminator locations respectively (expected for functional regulatory elements but not for random MITE insertions), the structural pre-existence of signals adjacent to proto-genes, and the compositional indistinguishability between boundary-proximal and dispersed runs (suggesting uniform distribution throughout the intergenic compartment rather than patchy MITE clustering). Whether MITEs seeded the AT regulatory landscape or whether they were themselves co-opted by a pre-existing regulatory system remains to be determined."

**RESOLUTION (draft4, 2026-03-21):** Fixed. MITE section in Discussion expanded into a full paragraph addressing: (1) strand-asymmetric positioning inconsistent with random MITE scars, (2) regulatory pre-existence adjacent to proto-genes argues systemic not insertion-driven, (3) compositional indistinguishability of boundary-proximal vs dispersed runs. Three alternative evolutionary scenarios for MITE-regulatory interaction proposed.

---

### 2.6 Monotone gradient Spearman ρ = -1.0 with N = 5 tiers — p-value not reported

**Where:** Results (proto-gene continuum section).

**The problem:** "Perfect monotonic gradient: GC decreases from 0.650 to 0.547 (ρ = -1.0)" is stated for the 5-tier trend. But with N = 5 data points (tiers, not genes), ρ = -1.0 has a specific one-tailed p-value. For N = 5, the number of possible orderings is 5! = 120. The probability of observing ρ = -1.0 (perfectly descending) by chance is 1/120 ≈ 0.0083. The probability of |ρ| = 1.0 (either perfectly ascending or descending) is 2/120 ≈ 0.017.

This is a marginally significant result with N = 5. Reporting ρ without the p-value leaves the reader unable to evaluate its significance. With two variables showing ρ = -1.0 (both GC and length), the evidence is stronger but the joint p-value should be computed (e.g., by Fisher's method combining the two p-values).

**Fix:** Add p-value explicitly: "GC content (Spearman ρ = -1.00, N = 5 tiers, exact p = 0.008 one-tailed) and gene length (ρ = -1.00, exact p = 0.008) decreased monotonically across confidence tiers." Then note that both metrics independently show the same gradient (Fisher's combined p ≈ 0.002), strengthening the inference.

**RESOLUTION (draft4, 2026-03-21):** Fixed. Table 2 caption updated with: "N=5, exact one-tailed p=0.0083 for each; Fisher combined p=0.00073 for the joint test." Note added that N=5 means each result is 1-in-120 by chance, and that biological significance lies in monotonicity across four independent measurements.

---

### 2.7 Independence assumption — within-genome gene comparisons

**Where:** Throughout all modules comparing gene subsets within P. salinus.

**The problem:** Virtually every statistical comparison in this paper treats the 1,430 P. salinus genes as 1,430 independent observations. They are not: all genes share a genome, meaning they share the same GC mutational background, the same host translational machinery, the same regulatory grammar, and the same evolutionary history. A Mann-Whitney test comparing GC3 between annotated and ORFan genes assumes the two gene sets are independent draws from populations — but both come from the same genome.

The practical consequence: reported p-values are likely anti-conservative (too small) because they ignore within-genome covariance. A codon usage comparison yielding p = 0.14 might be n.s. for both independence and non-independence reasons, but a comparison yielding p = 0.001 from n = 1,430 genes may be p = 0.05 after accounting for within-genome correlation.

This is a standard problem in comparative genomics that all computational biology papers face — you can't get around the fact that you're analyzing one genome. The correct response is to acknowledge it.

**Fix:** Add to Limitations: "All within-genome gene comparisons assume independence of individual gene statistics (GC content, CAI, RSCU distances, etc.). In reality, genes within a genome share evolutionary history, mutational background, and selective environment. The reported p-values for within-genome comparisons are therefore likely anti-conservative; they should be interpreted as measures of the consistency and magnitude of the differences observed rather than as calibrated false-positive rates. Key conclusions are based on large effect sizes (e.g., OR = 7.19 for the annotated-vs-hypothetical comparison) rather than on marginal p-values, providing some protection against this limitation."

---

### 2.8 Upstream 6-mer cosine similarity (0.9824): null distribution never reported

**Where:** Results (Regulatory signal section); Discussion.

**The problem:** "Annotated vs ORFan upstream 6-mer cosine similarity = 0.9824 — same regulatory grammar." But what would randomly selected intergenic sequences give? In a genome with 62% GC and extreme compositional uniformity, **random intergenic 6-mer profiles will also be highly similar** — the genome is so uniform that even random comparisons will give high cosine similarity.

If random intergenic comparisons give, say, 0.960 cosine similarity, then 0.9824 for annotated-vs-ORFan is only marginally above null and the "same regulatory grammar" conclusion is inflated. If random gives 0.85, then 0.9824 is genuinely high.

**Fix:** Report the null expectation: "The cosine similarity between random pairs of same-sized intergenic upstream sequences (n = 1,000 bootstrap replicates) was [X.XXXX ± SD], compared to the observed annotated-vs-ORFan similarity of 0.9824. The observed similarity is [above/at/below] the null distribution, indicating [regulated/compositionally driven] similarity in regulatory sequence content."

---

## TIER 3 — Polish and Tighten

These are quality-of-life improvements that strengthen reviewer confidence and paper precision. None are likely to generate major revisions on their own.

---

### 3.1 Bootstrapped CIs on major extrapolated estimates

The paper would be materially stronger with bootstrapped CIs on three estimates currently reported as point estimates:
- The extrapolated 91.4% (or whatever the correct) overall orphan rate
- The mean CAI difference (Pandoravirus vs Mimivirus): 0.130 difference needs a CI
- The "true intergenic drops from 32.5% to 23.8%": this is an approximation and should have an uncertainty range

These are straightforward to compute from existing data and would preempt reviewer requests.

---

### 3.2 Product-level hit rates — small-N disclaimer

The product-level gradient (ankyrin 42.4% → F-box 8.5% → Fascin 0%) is the most visually compelling finding in Figure 3. But some products have very small N: Ubiquitin n = 5 (40% hit rate = 2/5), BTB n = small. The figure caption should explicitly label products with n < 10 as "descriptive only."

More importantly: the denominator in ankyrin (56/132 = 42.4%) should be verified — does "132" represent all P. salinus genes with "ankyrin" in the product name? This is a large number; clarify how product assignment was made.

---

### 3.3 Hairpin excess (2.0pp, z = 3.6) — clarify "biologically small" with benchmarks

"Statistically significant but biologically small" for the 2.0pp hairpin excess is the correct interpretation but needs a benchmark. Biologically small compared to what? A brief reference: "The 2.0 percentage-point excess of palindromic hairpin sequences is considerably smaller than the 60–70% estimated from energy-based RNA folding by Legendre et al. (2018), suggesting that palindrome detection at the sequence level recovers a minor subset of the termination-compatible secondary structures identified by thermodynamic methods."

This also clarifies the reconciliation with Legendre more explicitly than the current framing.

---

### 3.4 Figure 4 Panel B — intergenic = 0% is definitional, not empirical

Figure 4B shows the "three-level maturation gradient: intergenic (0%) → hypothetical (3.0%) → annotated (18.2%)." The 0% for intergenic is **definitional** — intergenic regions are by definition not annotated genes, so they would not be queried in BLASTp. This point is in the "not genes" column of Table 4 but should be made explicit in the figure legend: "Intergenic regions were not queried by BLASTp; 0% represents the definitional baseline (non-genic sequence has no protein homologs), not an empirically measured absence."

Without this clarification, the gradient looks like three empirical data points when the first is actually a tautology.

---

### 3.5 Strand ratio r = -0.79 (N = 4) should be removed

**Where:** Module 2 discussion.

A Pearson correlation with N = 4 data points has df = 2. The critical r at p = 0.05 (two-tailed) with df = 2 is |r| ≥ 0.950. With r = -0.79, this correlation is **non-significant by its own standard** (p ≈ 0.21). Reporting it implies it is a meaningful result. The finding was already "downgraded from 'unusual' to 'needs more data'" but the r value should be removed entirely rather than left as a non-significant correlation that looks quantitative.

---

### 3.6 Phase 7b Proto-gene gradient: explicitly note the spurious proxy result

The paper (from memory) correctly uses the actual T4/T5 genes for the regulatory signal comparison (p = 0.14, p = 0.25), after discovering that the short-hypothetical proxy gave a spurious significant result (p = 0.0007). This correction is critical to the paper's integrity but should be explicitly described in the Methods as a methodological correction: "An initial analysis used short hypothetical genes (length < 500 bp) as a proxy for proto-gene candidates. This proxy showed significantly lower upstream AT content (p = 0.0007), but this result reflected the correlation between gene length and AT content rather than a genuine regulatory difference. Subsequent analysis used the actual T4/T5 proto-gene candidates, which showed no significant difference in upstream AT content (T4: p = 0.14, T5: p = 0.25), indicating that proto-genes already have regulatory signals indistinguishable from established genes."

This proactive disclosure of a corrected false positive is actually a strength — it shows analytical rigor. Don't hide it; feature it.

---

### 3.7 Hypothetical subsample: note random seed dependency

The 299 hypothetical genes queried used seeds 42 (first 100) and 43 (additional 200). The OR = 7.19 estimate is for this specific sample. A brief bootstrap or sensitivity note: "To assess sensitivity to the random subsample, we computed OR estimates for 100 additional random samples of n = 299 from the 902 hypothetical genes; all OR estimates fell within [X.X–X.X], confirming that the result is not sensitive to the specific subsample." If this has not been done, it should be — it takes 10 minutes of code and substantially strengthens the robustness claim.

---

### 3.8 ESM-2 GC confound: clarify which findings survive

The paper correctly notes "41% of PC1 correlates with GC content" as a limitation. But it doesn't explicitly tell the reader which structural findings are GC-confound-resistant and which are not. Add a brief summary:

Findings likely **robust** to GC confound:
- Proto-gene structural non-separability (p = 0.90): comparison is within the same embedding space, GC confound applies equally to both groups
- 23 ORFan singletons with structural similarity to known families (Fascin, ankyrin, kinase, Ubiquitin): these span the GC range

Findings potentially **confounded**:
- Low-GC outlier cluster (Finding 1): cluster may be defined by GC, not structure
- Proto-genes clustering with low-GC outliers (Finding 2): may follow from shared GC rather than shared structural theme

The Limitations section currently says "41% of PC1 variance correlates with GC content" but doesn't connect this to specific findings. Connect it.

---

### 3.9 The "p < 10^-6" in the abstract — give the exact value

The binomial test for 74% correct strand-aware positioning: with 3,559 boundary-associated runs (2,332 promoter + 1,227 terminator), p_observed = 0.74, p_null = 0.50, the exact binomial p-value is computable and should be in the range of 10^-100 or smaller. Report it exactly (or as < 10^-100 rather than < 10^-6). The 10^-6 threshold appears to be an approximation — the actual value is likely far more extreme.

---

### 3.10 Supplementary table listing all statistical tests

The statistics_audit.json already exists (140+ entries). Convert the 40-50 tests reported in the main text to a formatted Supplementary Table (test name, comparison, test type, test statistic, p-value, effect size). This serves three purposes:
1. Replaces the need for readers to hunt through the methods for each statistic
2. Demonstrates methodological transparency (GBE values this)
3. Makes the multiple testing situation visible and manageable

This could be Supplementary Table S10 (listed but currently described only as "statistics audit").

---

## Methodological Strengths (Explicitly Worth Noting in Cover Letter)

These should be highlighted in the cover letter or response to reviewers as evidence of methodological rigor:

1. **Multiple null calibrations** (Phases 7b and 8b) — most papers don't calibrate their methods against null distributions. The discovery that the background pairwise similarity (0.109) is a compositional artifact (null = 0.105, Cohen's d = 0.06) is methodologically exemplary.

2. **Taxonomy misclassification corrected** — catching and correcting the Micromonas/bacteria substring match bug before submission is exactly the kind of careful auditing that separates publishable from unreliable work.

3. **Connection failure / no-hit clarification** — correct denominator (521 not 528) reported and explained.

4. **E-value sensitivity analysis** (81.8% → 90.4% across four thresholds) — pre-empts the most common objection to BLASTp orphan rate estimates.

5. **Effect-size filtering in Module 3** — the obs/exp > 2.0 threshold alongside z > 3 is excellent practice for large-N genomic data. Acknowledging that "88.6% of significant k-mers lacked meaningful effect sizes" is unusually honest.

6. **Strand-aware boundary test** replacing the strand-naive proximity test (Phase 1b FIX 2) — catching and correcting this methodological flaw before the final analysis is the right scientific behavior.

7. **Fisher's exact with Haldane-corrected Woolf CI** — correct choice for the small-count 2×2 table (hypothetical vs annotated). Wilson score CI for proportions.

8. **Statistics audit JSON** (140+ entries) — this is exceptional for a single-investigator study. Cite it explicitly in the Data Availability statement and mention it in the cover letter.

9. **Prodigal as independent prediction tool** — using an ab initio predictor independently from the annotation being tested is genuine validation.

10. **Explicit CAI caveat** (Kazusa general table, not highly-expressed subset) — proactively acknowledging a methodological limitation before the reviewer catches it.

---

## Cross-Cutting Framework Issues (from Statistical Best Practices Framework)

These are the eight cross-cutting failures from the statistical best practices framework, evaluated for this paper:

| Failure | Status in this paper | Verdict |
|---------|---------------------|---------|
| **No uncertainty propagation** | CAI means, extrapolated orphan rates, effect sizes lack CIs | Partially addressed — needs bootstrapped CIs on 3 key estimates |
| **Annotation feedback loop** | Module 9 uses Legendre GFF (RNA-seq-informed, not raw reads); Module 6 uses NCBI annotation tiers | Acknowledged in limitations; cannot fully fix without raw TSS data |
| **Compositional data in Euclidean space** | GMM on GC/GC3/GC12/ENC/CAI without CLR transform | Minor issue — add to Limitations |
| **Pseudoreplication (within-genome)** | 1,430 genes treated as independent | Not acknowledged — add to Limitations (Tier 2.7 above) |
| **Causal overreach** | "Consistent with," "line of support" framing used throughout | Good — the careful language is correct and appropriate |
| **Effect size neglect** | Cohen's d reported for ESM-2; OR reported for BLAST; most comparisons lack effect sizes | Partially addressed — add to CAI comparison and proto-gene displacement |
| **Reference bias** | Only P. salinus analyzed; single-genome study | Appropriate — explicitly acknowledged as limitation |
| **Multiple testing without correction** | No FDR or family-wise correction | Needs explicit acknowledgment in Limitations (Tier 1.5 above) |

---

## Expected GBE Reviewer Questions — Pre-emptive Answers

These are the questions a GBE reviewer is most likely to ask. Having answers ready in the manuscript or in your response-to-reviewer template will accelerate the revision cycle.

**Q1: "The 92.5% ORFan rate in the abstract — how is this calculated? It doesn't match the reported BLASTp numbers."**
→ Fix before submission (Tier 1.1). Derive the number explicitly from the raw data.

**Q2: "You report p < 0.05 for many comparisons. With ~50 tests, how many false positives would you expect? Was any multiple testing correction applied?"**
→ Address in Limitations (Tier 1.5). Distinguish confirmatory from exploratory tests; emphasize convergent evidence from independent methods.

**Q3: "The neutrality plot slope is not reported. R² < 0.05 is insufficient to conclude selection is dominant."**
→ Fix in Results/Methods (Tier 1.2). Report slope and CI.

**Q4: "The comparison is based on 4 genomes. How do these findings generalize to other giant viruses or other Pandoravirus strains?"**
→ Pre-empt in Introduction and Limitations (Tier 2.1). Qualify all comparative superlatives.

**Q5: "You claim poly-A/T runs are regulatory signals, but MITE transposable elements are AT-rich and common in P. salinus intergenic space. Could MITEs explain the poly-A/T distribution?"**
→ Address in Discussion (Tier 2.5). Acknowledge the MITE alternative and explain why strand-asymmetry argues against pure MITE origin.

**Q6: "The proto-gene gradient could reflect annotation bias rather than biological maturation. Short, low-GC sequences are inherently harder to annotate."**
→ This is already partially addressed. Strengthen with the Phase 7c framing (Tier 2.3).

**Q7: "What is the CI on OR = 11.9 for the hypothetical spatial enrichment? With only 9 events, this estimate is highly uncertain."**
→ Report the CI before submission (Tier 1.4).

**Q8: "ESM-2 embeddings capture evolutionary/compositional patterns, not 3D structure. With 41% of PC1 driven by GC, how much of the structural signal is genuine?"**
→ Pre-empt with explicit statement of which findings are GC-confound-resistant (Tier 3.8).

**Q9: "The CAI comparison (0.787 vs 0.657) — is this difference statistically tested? What is the confidence interval?"**
→ Fix in Results (Tier 1.3).

**Q10: "You state T4 predictions are '10-15% of the way from random to genome-optimized.' What is the exact calculation?"**
→ Fix to match actual RSCU distance calculation (Tier 2.4). The actual displacement appears to be ~5-7% based on reported distances.

---

## Priority Action List (Ordered by Impact)

### Must do before submission:
1. Derive and verify the headline ORFan rate (currently 92.5%) from raw numbers; fix in abstract and throughout
2. Report neutrality plot slope + CI (not just R²)
3. Report Mann-Whitney p + bootstrap CI for CAI comparison (0.787 vs 0.657)
4. Compute and report Fisher's exact CI for OR = 11.9 (hypothetical spatial enrichment)
5. Add multiple testing acknowledgment to Limitations
6. Specify the exact test for proto-gene codon usage (T4/T5 vs random IG ORFs) in Methods
7. Specify the Prodigal overlap criterion precisely in Methods

### Should do before submission:
8. Qualify all comparative superlatives with "among the comparison genomes tested"
9. Add MITE-as-poly-A/T-source alternative to Discussion
10. Add GMM GC-correlation limitation note to Limitations
11. Add within-genome independence caveat to Limitations
12. Report the Spearman ρ = -1.0 p-value for N = 5 tiers (p ≈ 0.008, one-tailed)
13. Compute upstream 6-mer cosine similarity null (bootstrap from random intergenic sequences)
14. Verify and explicitly connect which ESM-2 findings are robust vs GC-confounded

### Polish if time allows:
15. Bootstrapped CIs on extrapolated estimates (orphan rate, CAI difference, intergenic fraction)
16. Product-level hit rates: label small-N products as "descriptive only"
17. Report exact binomial p-value for 74% strand-aware positioning
18. Explicitly describe the short-hyp proxy correction as a methodological correction in Methods
19. Bootstrap the hypothetical subsample OR estimate (resample from 902 hypothetical genes)
20. Remove the non-significant r = -0.79 strand ratio correlation (N = 4, meaningless)
21. Convert statistics_audit.json to Supplementary Table S10

---

## Summary Assessment

This paper is methodologically above average for single-investigator computational genomics and makes a genuine scientific contribution. The two-project structure (1.0 establishes compositional evidence, 2.0 adds regulatory and homology validation) is well-designed. The convergent evidence approach (12 independent analytical modules pointing at the same model) provides real protection against individual false positives.

The seven Tier 1 issues, if fixed, should make the paper defensible against the most likely major revision requests. The most important single fix is Item 1.1 (the 92.5% ORFan rate) — this is a headline number in the abstract that does not currently match the raw data, and a GBE reviewer will catch it within the first 15 minutes of reading.

**Recommended submission timeline:** Fix Tier 1 items (1–2 days), address high-priority Tier 2 items (2–3 days), then submit. GBE turnaround for initial review is typically 4–8 weeks.
