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
- **Relevance:**
  - Connects to the AT-content analysis — both studies reveal non-random compositional patterns
  - If an active DNA editing system avoids AGCT, it may also shape the poly-A/T run distribution
  - Supports the broader picture: Pandoravirus genomes have unusual sequence composition constraints beyond what random mutation/drift would produce

### 6. Sun et al. 2015 — "MITE transposable elements in the P. salinus genome"
- **Journal:** BMC Biology
- **Key findings:**
  - Miniature inverted-repeat transposable elements (MITEs) have colonized the P. salinus genome
  - MITEs reside in intergenic space — the same regions targeted by the compositional analysis
  - Represents mobile genetic element activity in the intergenic regulatory landscape
- **Relevance:**
  - MITEs occupy the intergenic space where poly-A/T regulatory signals reside
  - Could interact with poly-A/T signals: MITEs may disrupt, co-opt, or be recruited alongside AT-rich regulatory elements
  - Raises the question: are some poly-A/T runs MITE-derived, or do MITEs preferentially insert near existing AT-rich regions?

