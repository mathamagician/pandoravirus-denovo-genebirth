Dear Editors,

I am writing to submit the manuscript titled "Computational evidence for a distributed AT-based regulatory mechanism enabling de novo gene birth in *Pandoravirus salinus*" for consideration as a Research Article in *Genome Biology and Evolution*.

*Pandoravirus salinus* possesses the largest known viral genome (2.47 Mb), with 92.5% of its proteins lacking detectable homologs outside the Pandoraviridae (Philippe et al. 2013). Comparative genomic studies by the Claverie group have established that these genes likely arise de novo from non-coding sequence, but the genomic infrastructure enabling this process has remained unclear. Using a 12-module computational pipeline applied primarily to a single genome and validated across four additional Nucleocytoviricota, we identify a **distributed AT-based regulatory system** that pre-loads intergenic regions with promoter and terminator signals, resolving how proto-genes can be transcribed regardless of their genomic position. The strand-aware binomial test confirms 74% correct regulatory orientation in *P. salinus* (3,559/4,810; p ≈ 8 × 10⁻²⁵³) with same-genus replication in *P. dulcis* (61.8%; p ≈ 6 × 10⁻⁸⁰). Three secondary contributions support the central finding: a quantitative five-tier proto-gene continuum, an asymmetric four-stage maturation ordering (regulatory → structural → codon → GC), and cross-NCLDV validation showing singleton-dominated proteomes (64–91% under adaptive thresholds) and translational selection are broadly conserved across four families and 28–64% GC content.

I believe this manuscript is well suited for *GBE* because it addresses a fundamental question in genome evolution — how new genes arise from non-coding DNA — through computational methodology that complements the comparative-genomics approach. To our knowledge it would be the first NCLDV de novo gene birth study published in *GBE*, and it bridges two literatures (giant-virus genome evolution and de novo gene birth) that have developed largely independently. The work should be of interest to *GBE*'s readership in evolutionary genomics, giant virus biology, and de novo gene emergence.

We have made unusual efforts toward methodological transparency. The repository includes a per-claim statistics audit (`statistics_audit.json`) mapping every manuscript-cited number to its source computation and data file, all-vs-all BLASTp clustering of the 478-gene compositional family to validate the "distant paralogy" inference using sequence-level methods (Camacho et al. 2009), and full reproducibility of every figure from the deposited results. All analysis scripts, primary data, intermediate results, and the statistics audit file are publicly available at https://github.com/mathamagician/pandoravirus-denovo-genebirth and archived at https://doi.org/10.5281/zenodo.19347286.

This manuscript has not been published and is not under consideration elsewhere.

**Suggested reviewers** (we have no co-authorship or competing-interest relationships with any of the following):

1. **Frank O. Aylward** (Virginia Tech) — NCLDV taxonomy, ecology, and giant virus genome evolution. Author of the Aylward et al. 2021 *PLOS Biology* taxonomic framework cited in the manuscript.
2. **Aoife McLysaght** (Trinity College Dublin) — De novo gene birth methodology, ORFan analysis, and proto-gene continuum frameworks.
3. **Anne-Ruxandra Carvunis** (University of Pittsburgh) — De novo gene birth in yeast and across kingdoms; co-author of the canonical Van Oss & Carvunis 2019 *PLOS Genetics* review on de novo gene birth.

**Disclosure of AI tool use:** In accordance with *GBE* policy, I disclose that AI tools (Anthropic's Claude) were used to assist with code development, data-processing workflow design, analytical planning, and editorial refinement of manuscript text. This is also disclosed in the Materials and Methods section. All analyses, statistical methods, interpretations, and final manuscript content were reviewed and validated by the author, who takes full responsibility for the accuracy and integrity of the work.

Thank you for your consideration.

Sincerely,
Eddie Bradford
Bradford Genomics, San Diego, CA, USA
