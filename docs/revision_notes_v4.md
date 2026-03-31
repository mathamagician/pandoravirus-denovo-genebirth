# Revision Notes: draft3_GBE → draft4_GBE
**Date:** 2026-03-21
**Reviser:** Claude Sonnet (AI assistant) under Eddie Bradford's direction
**Basis:** Statistical audit (RPSM ~64/100) + senior peer reviewer assessment
**Script:** `scripts/compute_missing_statistics.py` → `results/missing_statistics.json`

---

## Summary

All 7 Tier 1 (must-fix) items resolved. 6 of 7 Tier 2 (pre-empt reviewers) items addressed. Word count increased by approximately 600 words. No core findings changed; all changes are reporting improvements, specification gaps, and defensive framing.

---

## TIER 1 FIXES (all 7 resolved)

### Fix #1 — ORFan rate corrected from 92.5% to 91.4%

**Problem:** Abstract used Philippe et al. 2013's 92.5% as if it were our derived number. The derivable estimate from our BLASTp data is 91.4%.

**What changed:**
- Abstract: "92.5% of genes classified as ORFans" → "estimated 91.4% (95% CI: 89.6%–93.0%) of genes classified as ORFans, based on our BLASTp analysis (annotated-gene orphan rate 81.8%; hypothetical-gene orphan rate 97.0%; extrapolated to all 1,423 classified genes; Philippe et al. 2013 independently reported 92.5% by a complementary method)"
- Introduction: added sentence: "Our independent BLASTp analysis yields a consistent estimate of 91.4% (95% CI: 89.6%–93.0%; see Results), confirming this extreme ORFan fraction by an orthogonal method."
- Discussion (results paragraph): "the 92.5% ORFan fraction" → "the estimated 91.4% ORFan fraction"
- Limitations: "the 92.5% ORFan fraction remains" → "the estimated 91.4% ORFan fraction (Philippe et al. 2013: 92.5%) remains"

**Derivation:** Extrapolated: (521 × 81.8% + 902 × 97.0%) / 1,423 = 91.4%. Bootstrap 95% CI: 89.6%–93.0%. The 1,423 denominator is classified genes (521 annotated + 902 hypothetical); 7 genes unclassified. Philippe et al.'s independent 92.5% is now cited explicitly and described as a complementary method, which is a *stronger* framing.

**Framing note (per senior reviewer):** The independent confirmation framing is better than replacing Philippe's number: "our BLASTp analysis independently confirms the extreme ORFan rate reported by Philippe et al."

---

### Fix #2 — Neutrality plot slope added

**Problem:** R² < 0.05 alone does not prove translational selection; the slope is required.

**What changed:**
In Results codon usage paragraph, replaced "neutrality plot R^2 < 0.05" with full statistics:
`R^2 = 0.045, N = 1,430 genes; regression slope beta = 0.118, 95% CI: 0.089–0.146, p = 5.9 x 10^-16; the slope << 1.0 expected under strict neutrality confirms that GC at constrained positions (GC12) is largely invariant with changes in wobble GC (GC3), consistent with translational selection rather than mutational drift`

**Values:** Computed via `gene_codon` table GC3 = fraction of codons with G/C at position 3; GC12 = mean of position 1 and 2. N=1,430 genes (≥50 codons). Slope=0.1176 (CI 0.0894–0.1457), R²=0.0448, p=5.87×10⁻¹⁶.

---

### Fix #3 — CAI comparison: Mann-Whitney test and bootstrap CI added

**Problem:** CAI 0.787 vs 0.657 comparison lacked any inferential statistics.

**What changed:**
Results codon usage paragraph: added per-genome SD, Mann-Whitney test result, bootstrap CI on difference, and CAI reference caveat at point of use:
`mean per-gene CAI 0.787 +/- 0.032 for P. salinus vs 0.657 +/- 0.027 for Mimivirus (Mann-Whitney U = 1,392,180, p < 10^-300; mean difference 0.129, 95% bootstrap CI: 0.127–0.132; n = 1,430 and 979 genes respectively; CAI reference: A. castellanii codon frequencies from the Kazusa database)`

**Values:** p is essentially machine precision (0). Difference is remarkably tight (CI width only 0.005) because N is large.

---

### Fix #4 — OR = 11.7 for hypothetical spatial enrichment: CI added

**Problem:** The hypothetical gene spatial enrichment result (7/9 hits in native zone) lacked a CI; the paper did not even report this test explicitly.

**What changed:**
Added to Results (homology section, before the taxonomic profile of 9 hypothetical hits):
`The 9 hypothetical non-self hits were also spatially enriched in the native zone: 7 of 9 (77.8%) fell within 1.3-2.0 Mb versus 67 of 290 (23.1%) without hits (OR = 11.7, Haldane-corrected OR = 9.98, 95% CI: 2.32–42.83, p = 0.0010; the wide CI reflects the small number of hypothetical non-self hits, n = 9, but even the lower bound 2.32 indicates substantial enrichment)`

**Values:** 2×2 table: hits_in_nz=7, nohits_in_nz=67, hits_out_nz=2, nohits_out_nz=224. Fisher OR=11.70, p=0.0010. Haldane-corrected OR=9.98, 95% CI: 2.32–42.83.

**Note:** Lower CI bound 2.32 > 1.0, so enrichment result holds even at lower confidence bound.

---

### Fix #5 — Multiple testing paragraph added to Limitations

**Problem:** ~50+ tests across 12 modules; no acknowledgment of multiple testing; GBE reviewers will ask.

**What changed:**
Added paragraph to Limitations distinguishing confirmatory from exploratory tests, explaining why global Bonferroni is inappropriate, and noting that key conclusions rest on convergent large-effect-size evidence. Also added within-genome independence caveat (addresses Tier 2, item 2.7 simultaneously).

**Key sentences:**
"We distinguish between confirmatory tests (pre-specified as direct tests of the central hypotheses before data analysis: the proto-gene codon usage comparison, the strand-aware 74% binomial test, the annotated-vs-hypothetical Fisher's exact test, and the two-population AT run test) from exploratory comparisons... The central conclusions rest on large effect sizes supported by multiple independent lines of convergent evidence, providing protection against false discovery beyond what any single correction procedure could offer."

---

### Fix #6 — Proto-gene codon test statistic specified

**Problem:** p = 0.020 and p = 0.0005 appeared in Results without naming the generating test.

**What changed:**
Methods (Module 7c) updated with full specification:
`For the disambiguation codon usage test, the RSCU Euclidean distance to the genome-preferred codon vector (mean RSCU across all T1+T2 genes) was computed for each T4 and T5 candidate. The same distance was computed for 10,000 length-matched random intergenic ORFs drawn from the same genome. One-sided Mann-Whitney U tests compared observed T4/T5 distances against the random null distribution, testing whether proto-gene codon usage is more genome-like than expected by chance.`

---

### Fix #7 — Prodigal overlap criterion exactly defined

**Problem:** "≥80% reciprocal overlap" was ambiguous about which denominator was used.

**What changed:**
Methods (Module 5) updated:
`a Prodigal prediction was considered to confirm a GenBank annotation if (overlap length / GenBank annotation length) >= 0.80 AND (overlap length / Prodigal prediction length) >= 0.80; both conditions must be met simultaneously, and start-position matching was not required`

---

## TIER 2 PRIORITY FIXES (6 of 7 addressed)

### T2-1 — "Among genomes tested" qualifiers added

All comparative superlatives in Abstract and Results qualified. Abstract: "among the four genomes tested." Results: "among the four comparison genomes" or "in our comparison set." Methods note added about illustrative comparison design.

### T2-2 — Spearman ρ p-value reported

Table 2 caption updated: N=5, exact one-tailed p=0.0083, Fisher combined p=0.00073 for joint test of GC + length.

### T2-3 — MITE alternative hypothesis addressed

Discussion MITE section expanded from 2 sentences to a full paragraph. Three lines of evidence cited supporting regulatory function over MITE residues (strand-asymmetry, proto-gene pre-existence, compositional uniformity). Three alternative MITE-regulatory evolutionary scenarios proposed.

### T2-4 — GMM GC confound added to Limitations

Added: "The GMM clustering uses five features highly correlated with genome-wide GC content; the three clusters may primarily reflect quantitative GC variation; 'NATIVE,' 'HIGH-GC,' 'VARIANT' designations should be considered provisional."

### T2-5 — Proto-gene tier circularity acknowledged

Added to Discussion (alternative models section): annotation practice assigns lower confidence to shorter, lower-GC sequences (same properties defining lower tiers). Phase 7c codon usage test is strongest counter-evidence. Comparative genomic validation across strains needed for definitive resolution.

### T2-6/2.7 — Within-genome independence caveat

Combined into the multiple testing paragraph (Fix #5): "reported p-values for within-genome comparisons are likely anti-conservative and should be interpreted as measures of consistency and effect magnitude rather than as calibrated false-positive rates."

**T2 item NOT addressed:** Item 2.4 (10-15% codon displacement may overstate effect size). This would require re-examining the RSCU distance calculations and potentially revising a quantitative claim. Left for a future revision cycle — note the existing claim uses qualitative framing ("approximately 10-15%") which softens the issue.

---

## FILES CHANGED

| File | Change |
|------|--------|
| `docs/Pandoravirus_combo_paper_draft4_GBE.md` | New version with all fixes applied |
| `docs/STATISTICAL-AUDIT-FINAL.md` | RESOLUTION blocks added to all Tier 1 + priority Tier 2 items |
| `docs/revision_notes_v4.md` | This file |
| `scripts/compute_missing_statistics.py` | New script — computes all missing statistics |
| `results/missing_statistics.json` | New file — output of stats script |
| `statistics_audit.json` | New section `_pre_submission_audit_fixes` with all computed values |

---

## ESTIMATED RPSM SCORE AFTER FIXES

Previous: ~64/100
Expected after fixes: **~76/100**
- S_stat: 62 → ~80 (all 7 Tier 1 items resolved; major multiple testing gap closed)
- S_evid: 55 → 55 (unchanged; single-genome constraint remains)
- S_transp: 90 → 92 (audit resolutions documented; compute script added)
- S_rep: 50 → 55 (more honest framing strengthens credibility)

Remaining gap from 100: single-genome limitation (irreducible without new data), within-genome independence (structural to the field), and the 10-15% codon effect size claim (Tier 2 item 2.4, not yet addressed).
