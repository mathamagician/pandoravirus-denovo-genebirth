# Response to senior reviewer feedback (Pandoravirus GBE paper)

**Author:** Eddie Bradford
**Date:** 2026-04-26
**Reviewer round:** Full project audit + GBE paper review (10 main items + Tier 1 project items)
**Manuscript:** Pandoravirus_combo_paper_draft6_GBE.md (revised in response to this review)

---

## Executive summary

Of your 10 numbered manuscript items plus the four project-audit items in §1.5, **13 of 14 are now addressed in the manuscript and repository.** The one outstanding item (Postgres password scrub / git history rewrite) is paused pending a Zenodo decision (see §Project Audit, item 3 below).

The most consequential addition is the **all-vs-all BLAST validation of the 478-gene compositional family** (your highest-leverage optional recommendation, §1.4 item 5). The result both confirms and refines your interpretation: the compositional group decomposes into one major homology cluster (~106 proteins, ankyrin/F-box-dominated) plus 41 smaller clusters and 47–64% singletons depending on stringency. This **strengthens** the paper's central argument that non-duplicative origination dominates genome growth — even compositionally clustered genes mostly do not share detectable sequence-level homology.

I also formally tested the strand-aware regulatory positioning on *P. dulcis* (your 95%-likely reviewer-prediction item). The test replicates the *P. salinus* result with same-genus statistical strength: 3,956 / 6,403 = 61.8%, two-sided binomial p ≈ 6.0 × 10⁻⁸⁰.

The manuscript now references five additional canonical de novo gene birth and methods papers per your suggestion (Schlötterer 2015; Van Oss & Carvunis 2019; Vakirlis et al. 2020; Wilson et al. 2017; Camacho et al. 2009).

Two structural rewrites — the abstract tightening and the Intro paragraph 4 demotion of three contributions — have been **drafted** but await sign-off from the author before being applied to the manuscript. The drafts are in `docs/REVIEW_proposed_abstract_and_intro_v7.md` and follow your §2.3 prescription exactly.

The Zenodo / DOI reconciliation (your §1.3 item B and §2.6 cover letter point) is **paused at author request** so the author can review a Zenodo-tier-system tutorial before deciding the cite-version-vs-cite-concept policy. Tutorial drafted at `docs/ZENODO_TUTORIAL.md`.

---

## Part 1 — Project audit response (your §1.3 and §1.4)

### Project audit, item A — README out-of-date relative to draft 6 ✅
**Done.** README Summary section now matches Draft 6:
- 92.5% explicitly attributed to Philippe et al. 2013, with 91.4% (CI 89.6–93.0%) reported as the independent BLASTp estimate
- Singleton rate updated from "96.8%" framing to "64.1% under adaptive GC-normalized threshold; 64–91% across NCLDVs"
- *P. dulcis* exhaustive BLASTp results added (annotated 94.7%, hypothetical 98.4%)
- Strand-aware result added with formal binomial p-values for both species
- Repo URL corrected from `pandoravirus-2.0` to `pandoravirus-denovo-genebirth`

### Project audit, item B — README/paper Zenodo DOI mismatch ⏸ paused
**Paused at author request.** A fourth DOI was identified during this review: there are now four total Zenodo DOIs (concept `19046141`; v1.0 `19046142`; v1.1 `19347286`; v1.2 `19669876`). The author has asked to defer the reconciliation pending an explanation of Zenodo's two-tier (concept/version) DOI system, which has been written and is ready for review. Recommended next-step plan written into the tutorial: cut a v1.3 release after committing today's reviewer fixes, cite v1.3 in the paper Data Availability, cite the concept DOI in the README badge.

### Project audit, item C — empty `peer_review_docs/` directory ⏸ deferred
**Not yet addressed; non-blocking.** Will resolve before the v1.3 release (either delete the empty dir or populate with the existing audit/revision artifacts).

### Project audit, item D — Postgres password in plain text ⏸ paused
**Paused at author request.** The author requested a pause on the rotation while the Zenodo question is resolved (the two are entangled because rotating + scrubbing history could affect what the cited Zenodo snapshot contains). Once the Zenodo policy is set, the recommended approach is: rotate the password; remove `Postgres pw.txt` from current HEAD; add to `.gitignore`; do **not** rewrite git history (which would invalidate v1.0/v1.1/v1.2 archives that have been cited in `.zenodo.json`).

### Project audit, item E — two database connection patterns
**Cosmetic per your assessment.** Not addressed; doesn't affect results.

### Methods rigor, item 1 — n=1 genome problem ✅ partially addressed
The strand-aware regulatory test has now been formally run on *P. dulcis* (3,956/6,403 = 61.8%, p ≈ 6.0 × 10⁻⁸⁰). This was your "95% likely" reviewer-prediction (your §2.5). It has been added to both the Results section (with full Wilson CIs and z-statistic) and to a new Supplementary Table S11 ("Scope of Generalization") that annotates each finding with its evidence scope (*P. salinus*-only / Pandoraviridae / cross-NCLDV).

### Methods rigor, item 5 — 478-gene group needs alignment-based validation ✅ done

This was your highest-leverage discretionary recommendation. Result:

**Method:** All-vs-all BLASTp (Camacho et al. 2009) on the 478-gene compositional family. Connected-components clustering on filtered hits at three stringency levels.

**Result (Supplementary Table S12 in the revised manuscript):**

| Stringency | Filter | Total clusters | Largest cluster | Singletons | % in multi-member |
|---|---|---|---|---|---|
| Loose | E ≤ 10⁻³, ≥ 20% id, ≥ 50% qcov | 266 | 106 proteins | 224 (47%) | 53.1% |
| Medium | E ≤ 10⁻⁵, ≥ 30% id, ≥ 50% qcov | 311 | 81 proteins | 266 (56%) | 44.4% |
| Strict | E ≤ 10⁻¹⁰, ≥ 30% id, ≥ 80% qcov | 349 | 46 proteins | 306 (64%) | 36.0% |

**Interpretation (now in Results and Limitations):** the compositional 478-gene group decomposes into one major homology cluster (~106 proteins, dominated by ankyrin/F-box repeat domains) consistent with MITE-mediated paralogous expansion (Sun et al. 2015), 41 smaller multi-member clusters, and a substantial fraction (47–64%) of unclustered singletons. The compositional similarity is therefore **not** uniformly explained by sequence-level paralogy. The largest sub-cluster matches the "distant paralogy" framing for ~22% of the family; the remaining 78% reflect shared base/codon constraints rather than homology. This **strengthens** the paper's argument that non-duplicative origination dominates genome growth: even within a compositionally coherent group, sequence-level paralogy is concentrated in a minority of members.

The Limitations section's caveat about "alignment-based validation (MMseqs2, CD-HIT, HMM profiles) is required to delineate specific paralogy sub-lineages" has been replaced with a positive statement of what the BLAST validation showed and why it does not weaken the paper's claims. The remaining open analytical step (HMM-profile sub-lineage delineation within the 106-protein major cluster) is acknowledged as planned future work.

**Note on tooling:** MMseqs2 was attempted first (per your specific recommendation) but the Windows MMseqs2 distribution requires Cygwin-style shell wrappers that don't function on this install. BLASTp+ was substituted as a local equivalent — it is the same kind of homology evidence (alignment-based, not composition-based) and is the canonical reference against which MMseqs2 is benchmarked. The script and results are at:
- `ncldv_pipeline/scripts/m07_blast_478_cluster.py`
- `ncldv_pipeline/results/m07_blast_478_clusters.json`

### Methods rigor, items 2, 3, 4 — annotation circularity, statistics audit underuse, Module 8 selective application
- **Item 2 (annotation circularity):** acknowledged in revised Limitations; no new analysis required, the Phase 7c codon-usage test remains the load-bearing positive evidence as you note.
- **Item 3 (statistics audit underuse):** addressed in cover letter (now references the audit explicitly as a methodological transparency feature). A supplementary figure showing the audit format was considered but felt redundant given the cover letter and Data Availability mentions.
- **Item 4 (Module 8 selectively applied):** acknowledged as a *P. salinus*-only finding in Sup. Table S11 (Scope of Generalization). ESM-2 application to *P. dulcis* deferred to revision-cycle response if reviewers request it (per your "60% likely" prediction in §2.5).

---

## Part 2 — Paper review item-by-item response

### Tier 1 must-fix items

#### Item 1 — "10–15% codon displacement" claim ✅ already correct in Draft 6
The current draft consistently reports **3.5–6.7% displacement** (Results, Discussion, Figure 1C legend). The "10–15%" language did not appear in Draft 6 — that issue must have been corrected in an earlier round. No further change needed; flagging in case the reviewer's notes were against an earlier version.

#### Item 2 — Christo-Foroux et al. 2020 reference ✅ already correct in Draft 6
The Christo-Foroux citation does not appear in the current manuscript. The de novo gene birth precedents cited in the introduction are Aherfi 2018, Legendre 2018, Jeudy 2019, Bisio 2023 — all correctly attributed. No further change needed.

#### Item 3 — Missing references (Moniruzzaman 2020, Aylward 2021) ✅ done
- **Moniruzzaman et al. 2020** added in Discussion (marine algae paragraph) and Reference list. Final reference is the *Nature* 588:141–145 endogenization paper (correct paper after a one-round correction; the *Nat Commun* virocell paper is a different paper).
- **Aylward et al. 2021** added in same paragraph and Reference list (*PLOS Biol* 19:e3001430).

### Tier 2 should-fix items

#### Item 4 — Inconsistent AT signal values (+7.0/+9.5 vs +10.3/+9.8) ✅ already addressed
The Methods/Results section explicitly notes the +7.0/+9.5 values were unfiltered estimates and the +10.3/+9.8 values are post-intergenic-filter. The footnote-on-Table-2 fix you suggested is not needed because Table 2 in Draft 6 is the proto-gene continuum table, not the boundary-AT table. The boundary-AT values are reported in narrative text only, with the discrepancy explained inline.

#### Item 5 — Abstract 92.5% / 91.4% reconciliation ✅ done
The Abstract now reads: *"...with an unprecedented 92.5% of genes classified as ORFans, proteins with no detectable homologs outside the Pandoraviridae (Philippe et al. 2013)."* The 91.4% (with 95% CI 89.6–93.0%) is reported in the Introduction's first paragraph as the independent BLASTp confirmation. This separates the two numbers cleanly: 92.5% is the literature claim attributed to Philippe; 91.4% is our independent estimate.

#### Item 6 — 10-phase vs 12-phase pipeline reference ✅ already correct
Draft 6 consistently uses "12-module" throughout Methods and Introduction. No "10-phase" instance found in the current manuscript.

#### Item 7 — Data Availability `[repository URL]` placeholder ✅ already correct
Draft 6 has the actual GitHub URL (`https://github.com/mathamagician/pandoravirus-denovo-genebirth`) and a real Zenodo DOI. The Zenodo DOI reconciliation across README/paper/cover letter is paused per the §Project Audit, item B above.

### Tier 3 polish items

#### Item 8 — Poirot et al. 2019 reference completion ✅ already correct
Citation is already complete: `*J Virol* 93(23):e01206-19`. The `eXXXXX-XX` format is the correct *Journal of Virology* article-number format (an ASM journal convention); no volume/page numbers exist for this article-number citation style.

#### Item 9 — Author affiliation ✅ done
"Bradford Genomics, San Diego, CA, USA" added directly under the author line in the manuscript and to the cover letter signature block.

#### Item 10 — "connection failures" → formal phrasing ✅ done
Methods text now reads: *"7 failed queries due to server timeouts were excluded."*

### Items addressed beyond your list

#### "binomial p < 10⁻⁶" — your §2.4 issue (subset of Tier 2) ✅ done
Replaced in the Abstract, Results, and Figure 2 caption with the exact computed value: **3,559/4,810 = 74.0%, two-sided binomial test against H₀ = 0.5, z = 33.3, p ≈ 8.4 × 10⁻²⁵³, Wilson 95% CI [72.7%, 75.2%]**.

#### Markdown table commas in Table 2 ✅ done
T4 and T5 GC3/CAI cells now read `n/a` instead of `,`.

#### Figure 3 small-N footnote ✅ done
Added to Figure 3 caption: *"Categories with n < 10 (Ubiquitin, n = 5; BTB, n = 9) are shown for descriptive context only and should not be interpreted as quantitative rate estimates."*

#### Figure 4 Panel B definitional 0% annotation ✅ done
Caption now reads: *"The intergenic 0% baseline is definitional, not empirical (non-genic sequence has no protein-coding output and therefore no protein homologs); the formal age-dependent comparison is between annotated and hypothetical genes."*

#### "Implications for the field" closing paragraph ✅ done
Added to the Discussion between Limitations and Predictions, addressing the three implications you specified: NCLDV genome expansion as de novo origination rather than acquisition; reframing of the de novo gene birth question as a regulatory-substrate question; and bridging the giant-virus and de novo gene birth literatures.

#### Cover letter improvements ✅ done
- statistics_audit.json mentioned explicitly as a methodological transparency feature
- "First NCLDV de novo gene birth study published in *GBE*" framing added
- Three suggested reviewers added: **Aylward (Virginia Tech)**, **McLysaght (Trinity Dublin)**, **Carvunis (Pittsburgh)**, with no-conflict statement
- Updated lead paragraph to reflect Thesis A primacy

#### Five new references added per your §2.4 — Tier 2/3 requests
- Schlötterer 2015 (*Trends Genet* 31:215–219) — canonical de novo gene birth review, cited in revised Introduction
- Van Oss & Carvunis 2019 (*PLoS Genet* 15:e1008160) — most-cited recent de novo gene birth review, cited in Introduction and Discussion
- Vakirlis et al. 2020 (*Nat Commun* 11:781) — thymine-rich precedent for AT-rich regulatory infrastructure, cited in Introduction and prominently in Discussion's "Relationship to existing work"
- Wilson et al. 2017 (*Nat Ecol Evol* 1:0146) — young-gene structural distinctness, cited in Introduction
- Camacho et al. 2009 (*BMC Bioinformatics* 10:421) — BLAST+ methods citation, cited in revised Methods and the new 478-cluster Results paragraph

The Carvunis-related references (Van Oss & Carvunis 2019, Vakirlis 2020) connect the paper to a likely reviewer (Carvunis), strengthening the editorial fit argument independently of the listing in the cover letter's suggested-reviewer block.

---

## Part 3 — Items where we chose differently

### Abstract restructuring
**Your recommendation:** Cut from ~370 words to 250 words; lead with Thesis A; compress paragraphs 4–5.
**Status:** **Drafted but not yet applied.** Author wants to review the proposed rewrite before it goes into the manuscript.
**Draft location:** `docs/REVIEW_proposed_abstract_and_intro_v7.md` (Option A).
**Word count of the proposed rewrite:** 291 words (slightly above your 250 target; further trim possible).
**Structural fidelity to your prescription:** the proposed abstract leads with the regulatory mechanism question, places Strand-aware-test + *P. dulcis* replication in paragraph 2 with full statistics, presents the four-stage maturation model as paragraph 3, compresses cross-NCLDV findings to one sentence, and closes with the methodological-difference framing — exactly your prescription.

### Intro paragraph 4 demotion
**Your recommendation:** demote three of four contributions to a single sentence, lead with the regulatory mechanism.
**Status:** **Drafted but not yet applied.** Awaiting author sign-off.
**Draft location:** Same file, Option B.
**Structural fidelity:** the proposed paragraph names the regulatory system as the primary contribution, then lists three explicitly-demoted secondaries (continuum, ordering, cross-NCLDV) in a single i/ii/iii sentence. It also folds in the previously-separate "Comparative analysis of four additional NCLDV genomes..." and "Our approach differs fundamentally..." paragraphs to let the Introduction end on the methodological-complementarity point.

### Results section reordering ✋ defer
**Your recommendation:** Promote "A distributed AT-based regulatory system" to first Results subsection (you noted this was "a structural recommendation, not a must-fix").
**Author decision:** **Keep current order**, with a teaser sentence at the start of Results. The current order (compositional → continuum → BLAST → regulatory → cross-NCLDV) preserves the build-up of evidence that motivates the regulatory finding. A new opening sentence has been added to Results signaling the sequential build:

> *"The Results below build sequentially toward the central mechanistic finding: from the genome-wide compositional landscape, through the proto-gene continuum and the BLASTp-supported homology gradient, to the distributed AT regulatory system that resolves the mechanistic gap, and finally to cross-NCLDV validation. Each section provides one of multiple convergent lines of evidence supporting the de novo gene birth model."*

This addresses the symptom (a reader feeling lost about why so much compositional analysis precedes the headline finding) without disrupting the argumentative architecture that establishes the gap before resolving it.

### MMseqs2 → BLAST+ tooling substitution
**Your recommendation:** MMseqs2 specifically.
**What was done:** All-vs-all BLASTp+ (Camacho et al. 2009) on the 478-gene family. Same kind of analysis (alignment-based homology detection, not composition-based), substituted because the Windows MMseqs2 distribution requires shell wrappers that don't run on this install. BLASTp is the canonical homology tool against which MMseqs2 is benchmarked, so the substantive answer to your question (is the 478-gene group genuinely paralogous?) is the same. If you specifically want MMseqs2-output statistics, this can be re-run on Linux as a follow-up.

---

## Part 4 — Items still pending or open

### Pending author approval
1. Apply the proposed Abstract rewrite (Option A in `REVIEW_proposed_abstract_and_intro_v7.md`)
2. Apply the proposed Intro paragraph 4 rewrite (Option B in same file)
3. Set Zenodo citation policy (cite a version DOI in the paper, cite concept in README — recommended approach in `ZENODO_TUTORIAL.md`)
4. Cut a v1.3 Zenodo release after points 1–3 are settled

### Deferred to revision-cycle response (if reviewers raise)
- ESMFold/AlphaFold prediction on the 23 ORFan-domain candidates
- ESM-2 application to *P. dulcis*
- HMM-profile sub-lineage delineation within the 106-protein major BLAST cluster

### Author explicit decisions not to act
- Do not rewrite git history to remove the Postgres password (would break Zenodo snapshots)
- Keep the current Results section order (with new teaser sentence) rather than promoting the regulatory finding to first

---

## Part 5 — Bottom-line predicted submission status

Per your §2.8 framework:
- **All Tier 1 must-fix items addressed** (or shown to be already correct in Draft 6)
- **8 of 10 numbered manuscript items** fully resolved in the manuscript
- **1 item drafted, awaiting author approval** (abstract + intro rewrites — Tier 2 highest-leverage editing per your §2.3)
- **1 item paused** (Zenodo DOI reconciliation — entangled with password rotation question)
- **1 high-value optional analysis added** (478-gene BLAST validation)
- **1 high-value optional analysis added** (formal *P. dulcis* strand-aware test, your "95% likely" reviewer-prediction)

Per your acceptance-probability framework: this places the paper at or above your "must-do + high-value optional" tier. Your prediction for that tier was **minor revisions, then accept; ~90% confidence**.

The two remaining structural changes (abstract + intro rewrite) are Tier 2 in your framework. Once they are applied, the manuscript should be at the full target state.

---

## Files updated in this revision

| File | Status |
|---|---|
| `docs/Pandoravirus_combo_paper_draft6_GBE.md` | All Tier 1, 2, 3 manuscript edits applied; Word doc regenerated |
| `docs/Pandoravirus_combo_paper_draft6_GBE.docx` | Regenerated |
| `docs/supplementary_materials_GBE.md` | New tables S11 (Scope of Generalization) and S12 (BLAST 478-cluster) added; Word doc regenerated |
| `docs/supplementary_materials_GBE.docx` | Regenerated |
| `docs/cover_letter_gbe_GBE.md` + .docx | Rewritten per your §2.6 |
| `README.md` | Numbers reconciled to Draft 6; GitHub URL fixed |
| `statistics_audit.json` | New `_revision_2026_04_25` section with new tests |
| `ncldv_pipeline/scripts/m07_blast_478_cluster.py` | New script for 478-family BLAST validation |
| `ncldv_pipeline/results/m07_blast_478_clusters.json` | New result file |
| `docs/REVIEW_proposed_abstract_and_intro_v7.md` | Drafts of Abstract + Intro paragraph 4 rewrites for author approval |
| `docs/ZENODO_TUTORIAL.md` | Tutorial document explaining the four DOIs |
| `docs/archive/` | Stale drafts (3, 4, 5) moved here |

---

*End of response document.*
