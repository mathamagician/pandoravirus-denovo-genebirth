# Literature Review — Pandoravirus Transcriptomics & De Novo Gene Birth

## Key Papers

### 1. Legendre et al. 2018 — "Diversity and evolution of the emerging Pandoraviridae family"
- **Journal:** Nature Communications 9:2285
- **DOI:** 10.1038/s41467-018-04698-4
- **PMC:** PMC5995976
- **Key findings:**
  - First comparative genomics + transcriptomics of Pandoraviridae
  - Strand-specific RNA-seq on P. salinus, P. dulcis, P. quercus, P. neocaledonia
  - mRNA extracted from synchronized infections, polyA+ enriched
  - 82.7-87% of pandoravirus genomes is transcribed
  - 70% of 3' transcript ends have palindromic hairpin termination signal
  - Median 3' UTR = 16 nt (very short)
  - Spliceosomal introns in 7.5-13% of genes
  - 157-268 long non-coding RNAs per strain
  - Open pangenome — each new strain adds ~60 new proteins
  - De novo gene creation proposed as main diversifying mechanism
  - Genome browser: http://www.igs.cnrs-mrs.fr/pandoraviruses/ (GBrowse)
  - Genomic accessions: KC977571 (P. salinus), KC977570 (P. dulcis), MG011689 (P. quercus), MG011690 (P. neocaledonia), MG011691 (P. macleodensis)
  - Proteomics: PXD008167 (PRIDE)
  - **RNA-seq raw data: NOT deposited in SRA/GEO/ENA** (only 17 DNA-seq entries exist for Pandoravirus in SRA)
  - Processed data available via GBrowse (GFF annotation + RNA-seq tracks)

### 2. Jeudy et al. 2019 — "Pandoravirus celtis illustrates the microevolution processes"
- **Journal:** Frontiers in Microbiology 10:430
- **PMC:** PMC6418002
- **Key findings:**
  - P. celtis and P. quercus: 96.7% identical genomes
  - ~30 genes unique to one or the other = novel genes post-divergence
  - All unique genes are strict ORFans with GC closer to intergenic regions
  - **Smoking-gun evidence:** pclt_cds_350 and pclt_cds_1084 are genes in P. celtis that correspond to intergenic regions in P. quercus
  - De novo gene creation from intergenic space confirmed
  - Slight negative selection on new genes
  - Novel proteins gain transcription from intergenic ORFs

### 3. Koonin/Abergel et al. 2023 — "Evolution of giant pandoravirus revealed by CRISPR/Cas9"
- **Journal:** Nature Communications 14:428
- **PMC:** PMC9879987
- **Key findings:**
  - CRISPR/Cas9 framework for Pandoravirus neocaledonia
  - Essential genes concentrated at 5' end of genome
  - ORFan/non-essential genes accumulate at 3' end
  - 25 copies of each chromosome; chain-reaction CRISPR to knockout all copies
  - Reduced core essential genome reminiscent of smaller ancestors
  - Genetic expansion increased genome robustness

### 4. Christo-Foroux et al. 2020 — "A large open pangenome and a small core genome"
- **Journal:** Frontiers in Microbiology 11:1486
- **Key findings:**
  - P. neocaledonia genes with orthologs enriched at 5' of genome
  - Strain-specific (ORFan) genes accumulate at 3' end
  - Biased expansion toward 3' end
  - Open pangenome with no convergence limit

### 5. Poirot et al. 2019 — "AGCT tetranucleotide anomaly in Pandoravirus genomes"
- **Journal:** Journal of Virology
- **Key findings:**
  - Complete absence of the AGCT 4-mer in some Pandoravirus strains
  - Points to a novel DNA editing or selection process operating on these genomes
  - Another example of unusual compositional signatures in Pandoravirus genomes
- **Relevance to our work:**
  - Connects to our AT-content analysis — both studies reveal non-random compositional patterns
  - If an active DNA editing system avoids AGCT, it may also shape the poly-A/T run distribution we observe
  - Supports the broader picture: Pandoravirus genomes have unusual sequence composition constraints beyond what random mutation/drift would produce

### 6. Sun et al. 2015 — "MITE transposable elements in the P. salinus genome"
- **Journal:** BMC Biology
- **Key findings:**
  - Miniature inverted-repeat transposable elements (MITEs) have colonized the P. salinus genome
  - MITEs reside in intergenic space — the same regions our analysis targets
  - Represents mobile genetic element activity in the intergenic regulatory landscape
- **Relevance to our work:**
  - MITEs are in the intergenic space where our poly-A/T regulatory signals reside
  - Could interact with poly-A/T signals: MITEs may disrupt, co-opt, or be recruited alongside AT-rich regulatory elements
  - Raises the question: are some poly-A/T runs MITE-derived, or do MITEs preferentially insert near existing AT-rich regions?
  - Worth noting in discussion as another force shaping the intergenic regulatory landscape

## Implications for Pandoravirus 2.0

### Phase 1 — Transcriptomic Validation
**Good news:** RNA-seq data exists and transcript boundaries have been mapped.
**Challenge:** Raw RNA-seq reads are NOT in SRA. Data is only accessible via:
1. GBrowse at IGS Marseille (GFF download + RNA-seq tracks)
2. Supplementary data from the 2018 paper
3. Potentially contacting authors for raw data

**Strategy options (in order of preference):**
1. **Download GFF from GBrowse** — contains reannotated gene coordinates with transcript evidence
2. **Extract transcript boundaries from supplementary tables** in Legendre 2018
3. **Use the hairpin signal findings** from the paper to compare with our poly-A/T predictions
4. **Contact authors** if needed for raw BAM/BigWig files

**Key comparison to make:**
- Our prediction: poly-A/T runs mark gene boundaries (promoters/terminators)
- Their finding: palindromic hairpin structures mark 3' transcript ends
- Question: Are these the same regions, overlapping, or complementary signals?
- The hairpin signals are at 70% of 3' ends — what about the other 30%? Are those our poly-A/T regions?

### Phase 2 — BLASTp Analysis
**Well-supported by literature:**
- Existing studies focused on pangenome comparison within Pandoraviridae
- Basic BLASTp against nr has NOT been systematically published for all 528 annotated genes
- This would be genuinely new information about evolutionary origins
- Could reveal whether annotated genes cluster with eukaryotic, bacterial, or viral homologs

### Framing for Publication
**Critical realization:** De novo gene birth in Pandoravirus is NOT our novel contribution (Claverie group proposed it first). Our novel contributions are:
1. **Independent confirmation** via different methodology (compositional analysis vs. comparative genomics)
2. **Maturation order discovery**: regulatory → structural → codon optimization → GC composition
3. **Mechanism**: distributed AT regulatory system pre-loading intergenic space
4. **Proto-gene continuum**: quantitative 5-tier gradient
5. **Regulatory pre-existence**: T4/T5 proto-genes already have mature regulatory signals

This should be framed as "computational evidence for the mechanism of de novo gene birth"
rather than "discovery of de novo gene birth."
