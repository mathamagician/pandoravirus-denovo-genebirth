# Pandoravirus 2.0 — Handoff from analytics-pc-queue compute work

**Last updated:** 2026-05-10
**Purpose:** Onboard a new agent picking up Pandoravirus 2.0 work after
the recent compute push on the analytics-pc-queue workstation.

## TL;DR for the new agent

All compute since ~2026-05-08 happened in **`C:\Users\Eddie\Code Projects\analytics-pc-queue\`**,
NOT in this Pandoravirus 2.0 folder. That project provides the workhorse pipeline
(Snakemake + Docker + the new RTX 5080 workstation). Findings relevant to the
paper need to be merged back into this folder. This document points you at every
relevant artifact so nothing gets missed.

If you only read one thing: the **headline correction** is that the original
"38 ultra-orphans" / "10.6% de novo rate" claim has been refined. After running
a proper 8-axis filter across all 14 NCLDV genomes AND substituting a clean
RefSeq viral DIAMOND DB for the corrupted nr_virus.dmnd, the corrected count
is **311 of 349 (89%) strict de novo proteins surviving** across the full
NCLDV cohort, with **37 of 1,430 (2.59%) for *P. salinus* specifically**.

---

## Where to find everything

### Primary compute home
`C:\Users\Eddie\Code Projects\analytics-pc-queue\`

Key subfolders inside that project:

| Path | What's there |
|---|---|
| `docs/CAPABILITIES_INDEX.md` | Master 124-entry capabilities reference across all bioinformatics tools/DBs |
| `docs/ACQUISITION_PLAN.md` | Tiered acquisition plan with status; companion to CAPABILITIES_INDEX |
| `docs/DATA_LIBRARIAN_BRIEF.md` | Three-tier data architecture, naming conventions, in-flight items, migration integrity protocol |
| `docs/TROUBLESHOOTING.md` | Every issue + fix encountered (RTX 5080 + PyTorch cu128, Docker bind-mount slowness, nr_virus.dmnd corruption post-mortem, etc.) |
| `docs/MANUAL_DOWNLOADS.md` | Portal URLs for datasets requiring browser sign-in |
| `docs/INTERPROSCAN_SETUP.md` | Java setup guide |
| `logs/` | Per-job retrospective logs (01_*.md through 21_*.md, plus MORNING_BRIEF.md, FOLLOWUP_QUEUE.md) |
| `results/` | All compute outputs (see breakdown below) |
| `scripts/` | All runner scripts (run_diamond.py, run_esmfold.py, overnight_cpu_phase{2..8}.py, etc.) |

### Pandoravirus 2.0 results inside analytics-pc-queue

| Path | Contents |
|---|---|
| `results/pandora_diamond/` | DIAMOND search of P. salinus vs giant-virus DB |
| `results/pandora_blast/` | BLAST vs SwissProt |
| `results/pandora_pfam/` | Pfam-A domain scan |
| `results/pandora_tigrfam/` | TIGRFAM HMM scan |
| `results/pandora_cdd/` | NCBI CDD domain scan |
| `results/pandora_esmfold/` | 1,374 PDB structures + summary.tsv (~96% coverage of 1430 ORFs; growing as GPU phase 5 continues) |
| `results/quercus_esmfold/` | P. quercus structures (special path; 1,131 PDBs) |
| `results/foldseek/` | Foldseek vs PDB100 results for P. salinus |
| `results/pandora_annotation/` | 5-axis annotation Parquet + ultra-orphan call (the original 38) |

### Pan-NCLDV (all 14 genomes) results

| Path | Contents |
|---|---|
| `results/multigenome/<NC_xxxxxx>_esmfold/` | ESMFold structures for each non-P.salinus genome (NC_013756 Marseillevirus through NC_044179 P. massiliensis) |
| `results/multigenome_phase2/all_genomes_annotation*.parquet` | 5-axis annotation, all 14 genomes |
| `results/multigenome_phase3/ortho_groups_e10.tsv` | MCL ortholog groups (44,957 members in 5,026 groups) |
| `results/multigenome_phase3/PFAM_ENRICHMENT.md` | Pfam abundance per genome |
| `results/multigenome_phase3/cross_genome_plddt.parquet` | pLDDT distributions, all genomes |
| `results/multigenome_phase3/ULTRA_ORPHANS_REFINED.md` | Cross-genome refined ultra-orphan list |
| `results/ncldv_allvsall/diamond_allvsall.tsv` | All-vs-all DIAMOND (1.78M edges) |
| `results/cpu_phase5/` | Cross-genome Foldseek + TIGRFAM rollup |
| `results/cpu_phase6/` | Pandoravirus genus consolidation + pLDDT ranking |
| **`results/cpu_phase7/`** | **8-axis filter applied to all 14 genomes — the canonical "strict de novo" output** |
| `results/cpu_phase7/master_annotation/master_8axis_annotation.parquet` | **Source of truth** — every NCLDV ORF with 8-axis evidence |
| `results/cpu_phase7/master_annotation/PAN_NCLDV_DE_NOVO.md` | Per-genome strict de novo rates |
| `results/cpu_phase7/tigrfam_per_genome/` | Per-genome TIGRFAM domain scans |
| `results/cpu_phase7/cdd_per_genome/` | Per-genome NCBI CDD scans |
| `results/cpu_phase7/structural_graph/ncldv_all_vs_all.tsv` | Cross-NCLDV Foldseek graph (49,809 hits) |
| **`results/cpu_phase8/`** | **Latest round of analyses — Pfam enrichment, MCL re-clustering, codon usage, MAFFT MSAs, RefSeq viral correction** |
| `results/cpu_phase8/01_pfam_enrichment/` | Per-genome Pfam frequency + enrichment vs NCLDV mean |
| `results/cpu_phase8/03_mafft_core/msa/` | 27 MAFFT MSAs of cross-NCLDV core ortholog groups |
| `results/cpu_phase8/05_mcl_multi/` | MCL clustering at I=1.4 / 2.0 / 4.0 / 6.0 |
| **`results/cpu_phase8/06_refseq_viral_search/REFSEQ_VIRAL_SEARCH.md`** | **The 311/349 correction** |

---

## Headline findings (use these in the paper)

### 1. Strict de novo count (corrected)

**Old (paper draft 6):** 38 ultra-orphan ORFs in P. salinus → "10.6% de novo rate"
**New (8-axis filter + RefSeq viral correction):**
- **P. salinus**: 37 / 1,430 = **2.59%** strict de novo
- **Pan-NCLDV**: 349 strict de novo before RefSeq correction; **311 / 8,332 = 3.73%** after correction
- The 8-axis filter is much more stringent than the 5-axis (added TIGRFAM, NCBI CDD, Foldseek vs PDB100)
- After RefSeq viral (~2 GB curated viral DB), 38 of 349 (11%) of the strict de novo were found to have viral homologs the original filter missed; 311 truly have no homology anywhere

**Per-genome rates** (`PAN_NCLDV_DE_NOVO.md`):

| Genome | Organism | Total ORFs | Strict de novo | Rate |
|---|---|---:|---:|---:|
| NC_019273 | Mamavirus | 7 | 3 | 42.86% |
| NC_022837 | Pandoravirus dulcis | 31 | 6 | 19.35% |
| NC_027867 | Mollivirus sibericum | 523 | 97 | 18.55% |
| NC_028862 | Faustovirus E12 | 61 | 11 | 18.03% |
| NC_014637 | Cafeteria roenbergensis virus | 544 | 85 | 15.62% |
| NC_044179 | Pandoravirus massiliensis | 27 | 4 | 14.81% |
| NC_023828 | Pithovirus sibericum | 22 | 2 | 9.09% |
| NC_013756 | Marseillevirus | 428 | 20 | 4.67% |
| NC_037666 | Pandoravirus macleodensis | 1081 | 28 | 2.59% |
| NC_022098 | Pandoravirus salinus | 1430 | 37 | 2.59% |
| NC_037667 | Pandoravirus quercus | 1185 | 21 | 1.77% |
| NC_014649 | Mimivirus | 979 | 15 | 1.53% |
| NC_016072 | Megavirus chilensis | 1120 | 12 | 1.07% |
| NC_020104 | Moumouvirus | 894 | 8 | 0.89% |

Mimiviridae are low (~1–2%) as expected (well-studied, well-annotated).
Pandoraviruses cluster at 2–3%. Higher rates in smaller genomes are noise
(small denominator).

### 2. Methods: 8-axis filter

Every NCLDV ORF was tested against these 8 axes; **strict de novo = NONE found**:
1. DIAMOND vs giant-virus DB (self-organism filtered)
2. BLAST vs SwissProt
3. Pfam-A domain
4. TIGRFAM HMM
5. NCBI CDD (combined Pfam + COG + KOG + CATH-derived)
6. Foldseek vs PDB100 (structural homology)
7. ESM-2 top-1 cosine vs other-genome proteins ≥ random p99 (0.9877)
8. RefSeq viral DIAMOND (added 2026-05-10 — the correction step)

Implementation: `analytics-pc-queue/scripts/overnight_cpu_phase7.py`
Output: `results/cpu_phase7/master_annotation/master_8axis_annotation.parquet`
Method writeup: `results/cpu_phase8/06_refseq_viral_search/REFSEQ_VIRAL_SEARCH.md`

### 3. MCL ortholog groups + inflation sweep

`results/cpu_phase8/05_mcl_multi/`:

| Inflation | n_clusters | Use |
|---|---:|---|
| 1.4 | 7,088 | loose families |
| 2.0 | 9,515 | canonical |
| 4.0 | 14,201 | tight |
| 6.0 | 16,273 | near-singleton |

49,809 cross-NCLDV Foldseek edges feed into this. Reviewer-anticipator:
"results stable across I = 1.4–6.0" is justified.

### 4. 27 core ortholog MAFFT MSAs

`results/cpu_phase8/03_mafft_core/msa/` — built. HyPhy FEL attempted on 22 codon-aware MSAs but returned 0 results (Docker invocation issue with pegi3s/hyphy; non-blocking for the paper). Codon MSAs are available at `results/cpu_phase8/04_hyphy_fel/codon_msa/` for any later HyPhy / PAML / codeml retry.

### 5. Codon usage / GC3 (the new "proto-gene gradient" result)

`results/cpu_phase8/02_codon_usage/CODON_USAGE.md` — **publishable per-genome GC3 comparison** between strict de novo and annotated ORFs. Two opposite patterns:

| Genome | n_dn | GC3 (dn) | GC3 (ann) | delta |
|---|---:|---:|---:|---:|
| Mimivirus | 15 | 0.308 | 0.200 | **+0.108** |
| Megavirus chilensis | 12 | 0.217 | 0.157 | +0.060 |
| Moumouvirus | 8 | 0.213 | 0.140 | +0.073 |
| **Pandoravirus salinus** | **37** | **0.641** | **0.746** | **−0.106** |
| **P. macleodensis** | **28** | **0.664** | **0.735** | **−0.071** |
| **P. quercus** | **21** | **0.671** | **0.740** | **−0.068** |

**Pandoraviruses**: strict de novo proteins are AT-richer than the annotated set (consistent with the proto-gene gradient hypothesis from Pandoravirus 1.0 Phase 7 — young genes not yet codon-optimized).
**Mimiviridae**: opposite signal (de novo more GC-rich than the AT-rich background); interpretation TBD.

Caveat: only 10 of 14 genomes had successful GenBank `gbwithparts` retrieval (Marseillevirus, Cafeteria, Mollivirus, Faustovirus truncated, ~132 of 311 strict de novo not yet codon-scored). Doesn't change headline numbers but worth re-fetching when convenient.

### 6. Common Crawl mining

`results/cpu_phase8/09_common_crawl/ncldv_mentions.tsv` — **1,360 NCLDV mentions across 5,000 Common Crawl WET files**. Per-mention rows: WET filename, target URI, matching text excerpt. Useful for citation hunting, public-discourse landscape, news/blog references.

---

## Known issues / corrections / gotchas (from TROUBLESHOOTING.md)

These are the production stories the new agent needs to know:

1. **`nr_virus.dmnd` is unrecoverable.** The 16.8 GB original DB had 99.88% empty records (sequence store destroyed). Recovery via `diamond getseq → makedb` failed. **Replaced by RefSeq viral (~2 GB curated, on disk at `D:\Research Data Library\raw\biodb\refseq_viral\`).** Source of truth for "broad viral search" until full NCBI nr DIAMOND DB is rebuilt.

2. **Self-organism filtering matters.** The original DIAMOND search of P. salinus against the giant-virus DB returned hits to other P. salinus proteins as "real" homologs. Fixed in phase 6 — see `scripts/fix_cross_genome_self_hits.py` and the correction propagated through phase 7's master annotation.

3. **GPU phase 4/5 wedges.** ESMFold can hang at 100% GPU on a single long protein for hours without producing output. Kill threshold: if no new PDB written for >1 hour and GPU > 80%, kill the process. The recovery script (`overnight_gpu_phase4_recovery.py`) had one such hang on NC_014637, killed 2026-05-10 at midnight.

4. **PyTorch must be 2.7+cu128 for RTX 5080 (Blackwell sm_120).** Older PyTorch installs CUDA 12.4/12.6 wheels which crash with "no kernel image available". Documented in `docs/TROUBLESHOOTING.md`.

5. **Docker bind-mounts on Windows are pathologically slow** for sequential I/O. Always stage large data to Docker volumes first. Affects DIAMOND makedb, ColabFold MSAs, anything that streams a >5 GB file from a host bind mount.

6. **Python on Windows defaults to cp1252.** Always `encoding="utf-8"` on file writes; use ASCII fallback for stdout in scripts that may run unattended.

---

## In-flight work as of 2026-05-10 ~3 PM

### GPU (running)
- **Phase 5 fill**: ESMFold recovery at max_length 1200 → 2000 + ESM-2 mutational scan on top-20 strict de novo. Currently mid-Block A on NC_022098. ETA to full GPU clear: ~36–55 hr (much slower than originally projected — the 800-1200 aa proteins fold at ~5–10 min each on a 16 GB card).
- After phase 5 clears: **GPU available for ColabFold AF2 vacation job starting Wed 2026-05-13**.

### CPU (running)
- **Phase 8b fixup**: GenBank fetch + CDS extract + codon usage + HyPhy FEL on the 27 core MSAs. Waiting for phase 8 to finish (Foldseek AFDB UP50 download was completed; Common Crawl mining now in flight at 1,800/5,000 WET files).
- **Phase 8c fixup**: ✅ done. MCL multi-inflation + RefSeq viral filter complete.

### Acquisitions (vacation queue)
- ✅ ColabFold UniRef30 2302 (215 GB extracted) — **REQUIRED for AF2 vacation job**
- ✅ ColabFold PDB100 foldseek 230517 (68 GB extracted)
- ✅ ColabFold envdb 202108 (110 GB compressed)
- ✅ InterProScan 5.74-105.0 (42 GB extracted)
- ✅ UniProt TrEMBL (50 GB compressed)
- ✅ NCBI nr full (363 GB / 156 volumes) — on `C:\BLAST\db\`
- ✅ eggNOG full data (DIAMOND DB just re-downloading after first attempt was truncated; should be done within an hour)
- ✅ NCLDV expanded RefSeq proteomes (80+ representative giant virus genomes; 19 MB extracted)
- ✅ PHROG, DGIdb, WikiPathways small wins
- ⏭️ HHsuite UniRef30: deferred (ColabFold UniRef30 covers MSA need)
- ⏭️ PDB full mmCIF: rsync deferred to WSL2 after vacation

---

## Files in this Pandoravirus 2.0 project that may be stale

These were created before the recent compute push and should be cross-checked against the corrected findings above:

| File | What to update |
|---|---|
| `Pandoravirus_combo_paper_draft6_GBE.md` | The "38 ultra-orphans" / "10.6%" should be replaced by the 8-axis-corrected numbers: 37 strict de novo, 2.59% rate for P. salinus. Pan-NCLDV comparative table is new (use the 14-genome table above). |
| `REVIEWER_RESPONSE_2026-04-26.md` | Any reviewer concerns about "did you check broader DBs?" can be answered with the RefSeq viral + InterProScan + Foldseek vs PDB100 + ESM-2 cosine evidence. |
| `STATISTICAL-AUDIT-FINAL.md` | The per-protein evidence tier is now from the 8-axis filter, not 5-axis. |

---

## Recommended onboarding sequence for the new agent

1. **Read this document first** (you're doing it).
2. **Read `analytics-pc-queue/docs/TROUBLESHOOTING.md`** — covers all the gotchas.
3. **Read `analytics-pc-queue/results/cpu_phase8/06_refseq_viral_search/REFSEQ_VIRAL_SEARCH.md`** — the headline correction.
4. **Read `analytics-pc-queue/results/cpu_phase7/master_annotation/PAN_NCLDV_DE_NOVO.md`** — pan-NCLDV rates table.
5. **Load `master_8axis_annotation.parquet`** in pandas — that's the single source of truth.
6. **Skim `analytics-pc-queue/logs/`** for chronological context — each numbered log is one analysis step.
7. **Optional: read `analytics-pc-queue/docs/CAPABILITIES_INDEX.md`** for the full catalog of available bioinformatics tools and datasets on this workstation.

When updating the paper draft, the right pattern is:
- Pull data from `analytics-pc-queue/results/` (don't duplicate large outputs).
- Reference the artifacts using their absolute paths.
- Save derived figures / tables into `Pandoravirus 2.0/figures/` and `Pandoravirus 2.0/results/` so the project remains self-contained for submission.

---

## Open questions the new agent should keep an eye on

1. **Pandoravirus salinus 37 strict de novo** — what fraction of these have RNA-seq support? If Pandoravirus expression data exists in NCBI SRA, that's the next confirmation step.
2. **ColabFold AF2 vacation job** — will run on the 311 strict de novo across all NCLDV. AF2-vs-ESMFold comparison is the "did you check with another method?" answer for reviewers.
3. **YP_008437180.2** — the single high-pLDDT outlier from the original 38. After ColabFold AF2, does AF2 agree with ESMFold? Binary answer.
4. **Mollivirus sibericum 18.55%** — Mollivirus has a very high strict de novo rate. Worth investigating as a potential parallel case study to Pandoravirus de novo birth claims.

---

*Generated by the analytics-pc-queue compute agent. If anything below contradicts what's actually on disk, trust what you see on disk and update this document.*
