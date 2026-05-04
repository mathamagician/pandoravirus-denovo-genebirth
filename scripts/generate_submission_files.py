"""
Generate GBE submission files:
1. Manuscript .docx (user adds line numbers in Word, exports to PDF)
2. Cover letter .docx (user exports to PDF)
3. Supplementary materials .docx (tables S1-S10 + figures S1-S6, user exports to PDF)
"""

import pypandoc
import os
import csv
from pathlib import Path

BASE = Path(r"C:\Users\Eddie\Code Projects\Pandoravirus 2.0")
DOCS = BASE / "docs"
SUPP = BASE / "supplementary"
FIGS_2 = BASE / "figures"
FIGS_1 = Path(r"C:\Users\Eddie\Code Projects\Pandoravirus\figures")

# ============================================================
# 1. Manuscript .docx
# ============================================================
print("=== Generating manuscript .docx ===")
manuscript_md = str(DOCS / "Pandoravirus_combo_paper_draft3_GBE.md")
manuscript_docx = str(DOCS / "Pandoravirus_combo_paper_draft3_GBE.docx")

pypandoc.convert_file(
    manuscript_md, 'docx',
    outputfile=manuscript_docx,
    extra_args=['--wrap=none']
)
print(f"  -> {manuscript_docx}")

# ============================================================
# 2. Cover letter .docx
# ============================================================
print("=== Generating cover letter .docx ===")
cover_md = str(DOCS / "cover_letter_gbe_GBE.md")
cover_docx = str(DOCS / "cover_letter_gbe_GBE.docx")

pypandoc.convert_file(
    cover_md, 'docx',
    outputfile=cover_docx,
    extra_args=['--wrap=none']
)
print(f"  -> {cover_docx}")

# ============================================================
# 3. Supplementary materials — build markdown then convert
# ============================================================
print("=== Generating supplementary materials ===")

supp_md_parts = []
supp_md_parts.append("# Supplementary Materials\n\n")
supp_md_parts.append("**Computational evidence for a distributed AT-based regulatory mechanism enabling de novo gene birth in *Pandoravirus salinus***\n\n")
supp_md_parts.append("Eddie Bradford\n\n")
supp_md_parts.append("---\n\n")

# --- Supplementary Tables ---

table_info = [
    ("S1", "table_s1_codon_usage.csv",
     "Codon usage comparison between annotated and ORFan genes: RSCU, GC3, GC12, ENC, and CAI with statistical tests."),
    ("S2", "table_s2_structural_families.csv",
     "ESM-2 structural family assignments for 170 representative proteins, including tier, compositional cluster, and nearest-neighbor similarity."),
    ("S3", "table_s3_gene_families.csv",
     "Gene family detection results: all 14 compositional families with member counts, mean similarity, and spatial distribution."),
    ("S4", "table_s4_product_hit_rates.csv",
     "Product-level non-self BLASTp hit rates for all annotated protein categories."),
    ("S5", "table_s5_evalue_sensitivity.csv",
     "E-value sensitivity analysis for annotated gene non-self hit rates across four thresholds."),
    ("S6", "table_s6_nonself_hits_annotated.csv",
     "Complete taxonomy of all 95 non-self BLAST hits with organism, identity, coverage, and genomic position."),
    ("S7", "table_s7_nonself_hits_hypothetical.csv",
     "All 9 hypothetical gene non-self hits with taxonomy, identity, and native zone membership."),
    ("S8", "table_s8_crossgenome_cv.csv",
     "Cross-genome compositional homogeneity: full GC, GC3, GC12, ENC, and CAI coefficient of variation comparison."),
    ("S9", "table_s9_hypothesis_matrix.csv",
     "Hypothesis testing matrix: evidentiary status of four genome growth models (HGT, gene duplication, mobile element expansion, de novo gene birth) across all 12 modules."),
    ("S10", "table_s10_at_boundary_signals.csv",
     "Gene boundary AT signals across genomes: AT content, start and end delta-AT values, and poly-A/T intergenic fraction for Lambda, Mimivirus, and P. salinus."),
]

for label, filename, caption in table_info:
    filepath = SUPP / filename
    supp_md_parts.append(f"## Supplementary Table {label}\n\n")
    supp_md_parts.append(f"**Table {label}.** {caption}\n\n")

    # Read CSV and convert to markdown table
    with open(filepath, 'r', encoding='utf-8') as f:
        reader = csv.reader(f)
        rows = list(reader)

    if rows:
        # Header
        header = rows[0]
        supp_md_parts.append("| " + " | ".join(header) + " |\n")
        supp_md_parts.append("| " + " | ".join(["---"] * len(header)) + " |\n")
        # Data rows (limit to first 50 for very large tables like S1, S4, S6)
        data_rows = rows[1:]
        max_rows = 50 if label in ("S1", "S4", "S6") else len(data_rows)
        for row in data_rows[:max_rows]:
            # Escape pipe characters in cell values
            escaped = [cell.replace("|", "\\|") for cell in row]
            supp_md_parts.append("| " + " | ".join(escaped) + " |\n")
        if len(data_rows) > max_rows:
            supp_md_parts.append(f"\n*Table truncated; full data ({len(data_rows)} rows) available in the repository CSV file.*\n")

    supp_md_parts.append("\n---\n\n")

# --- Supplementary Figures ---

# Map supplementary figures to source files
supp_figures = [
    ("S1", "Genome architecture: sliding-window GC profiles, gene density maps, GC skew plots, and strand maps for all four genomes in the comparison ladder.",
     [
         (FIGS_1 / "phix174_genome_map.png", "PhiX174 genome map"),
         (FIGS_1 / "lambda_genome_map.png", "Lambda genome map"),
         (FIGS_1 / "mimivirus_genome_map.png", "Mimivirus genome map"),
         (FIGS_1 / "pandoravirus_genome_map.png", "P. salinus genome map"),
     ]),
    ("S2", "GMM clustering of P. salinus genes: three compositional clusters with native zone (1.3-2.0 Mb) overlay, showing cluster composition shift within the zone.",
     [(FIGS_1 / "pandoravirus_origin_classification.png", "GMM clustering")]),
    ("S3", "Pairwise protein 3-mer cosine similarity distributions: P. salinus vs Mimivirus, with composition-matched null model overlay.",
     [(FIGS_1 / "gene_family_similarity_distributions.png", "Similarity distributions"),
      (FIGS_1 / "null_similarity_calibration.png", "Null calibration")]),
    ("S4", "ESM-2 PCA embedding space showing structural family assignments, proto-gene positions, GC content overlay, and null calibration.",
     [(FIGS_1 / "phase8_structure_prediction.png", "Structure prediction"),
      (FIGS_1 / "phase8b_null_calibration.png", "Structure null calibration")]),
    ("S5", "Full 7-panel Module 11 BLASTp analysis figure.",
     [(FIGS_2 / "phase2_blastp_homology.png", "BLASTp analysis")]),
    ("S6", "Module 10 corrected analysis 6-panel figure.",
     [(FIGS_2 / "phase1b_corrected.png", "Transcriptomic validation")]),
]

for label, caption, image_files in supp_figures:
    supp_md_parts.append(f"## Supplementary Figure {label}\n\n")
    supp_md_parts.append(f"**Figure {label}.** {caption}\n\n")
    for img_path, alt in image_files:
        if img_path.exists():
            # Use forward slashes for pandoc compatibility
            img_str = str(img_path).replace("\\", "/")
            supp_md_parts.append(f"![{alt}]({img_str})\n\n")
        else:
            supp_md_parts.append(f"*[Image not found: {img_path.name}]*\n\n")
    supp_md_parts.append("---\n\n")

# Write combined supplementary markdown
supp_md_path = DOCS / "supplementary_materials_GBE.md"
with open(supp_md_path, 'w', encoding='utf-8') as f:
    f.write("".join(supp_md_parts))

# Convert to docx
supp_docx = str(DOCS / "supplementary_materials_GBE.docx")
pypandoc.convert_file(
    str(supp_md_path), 'docx',
    outputfile=supp_docx,
    extra_args=['--wrap=none']
)
print(f"  -> {supp_docx}")

print("\n=== All submission files generated ===")
print(f"  1. Manuscript:      {manuscript_docx}")
print(f"  2. Cover letter:    {cover_docx}")
print(f"  3. Supplementary:   {supp_docx}")
print("\nNext steps:")
print("  - Open manuscript .docx in Word, add line numbers (Layout > Line Numbers > Continuous)")
print("  - Fix any superscript formatting (10^-6 -> 10⁻⁶, etc.)")
print("  - Export all three as PDF")
print("  - Upload to https://academic.oup.com/gbe")
